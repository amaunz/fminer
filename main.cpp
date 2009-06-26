// fminer.cpp
// © 2008 by Andreas Maunz, andreas@maunz.de, jun 2008

/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <getopt.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <cassert>

#include "fminer.h"

using namespace std;
using namespace OpenBabel;

extern Statistics* statistics;

Fminer* fminer;

// helper routines
void puti ( FILE *f, int i ) {
  char array[100]; 
  int k = 0;
  do {
    array[k] = ( i % 10 ) + '0';
    i /= 10;
    k++; 
  }
  while ( i != 0 ); 
  do {
    k--;
    putc ( array[k], f );
  } while ( k );
}

void remove_dos_cr(string* str) {
	string nl = "\r"; 
	for (string::size_type i = str->find(nl); i!=string::npos; i=str->find(nl)) str->erase(i,1); // erase dos cr
}

void read_gsp (char* graph_file) {
    FILE *input = fopen (graph_file, "r");
    if (!input) {
        cerr << "Error opening file '" << graph_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }
    fminer->ReadGsp(input);
}

void read_smi (char* graph_file) {
    Tid tree_id = 0;

    ifstream input;
    string line;
    int line_nr=0;
    string tmp_field="";

    // Open input stream
    input.open(graph_file);
    if (!input) {
        cerr << "Error opening file '" << graph_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }

    // draw line from input
    int field_nr = 2; // initialize to 'valid' pattern
    while (getline(input, line)) {
        line_nr++;
        istringstream iss(line); // line as iss
        string id = "";

        // if previous line has had two entries, i.e. format was correct, then re-initialize
        if (field_nr==2) field_nr=0;
        else {
            cerr << "Error! Line no. " << line_nr-1 << " in file '" << graph_file << "' is not in correct input format. Please refer to README." << endl;
            exit(1);
        }

        while(getline(iss, tmp_field, '\t')) {  // split at tabs
            if (field_nr == 0) {        // ID
                tree_id = (unsigned int) atoi (tmp_field.c_str());
                if (tree_id == 0) { 
                    cerr << "Error! Invalid ID: '" << tmp_field << "' in file  '" << graph_file << "', line " << line_nr+1 << "." << endl; 
                    exit(1); 
                }
                field_nr=1;
            }
            else if (field_nr == 1) {   // SMILES
                fminer->AddCompound (tmp_field , tree_id);
                field_nr=2;
            }
			else {
                cerr << "Error! Line no. " << line_nr << " in file '" << graph_file << "' is not in correct input format. Please refer to README." << endl;
				exit(1);
			}
        }
    }
    cerr << fminer->GetNoCompounds() << " compounds" << endl;
    input.close();

}

void read_act (char* act_file) {
    ifstream input;
    string line;
    string tmp_field;
    string act_name;
    Tid tid=0;
    unsigned int line_nr = 0;
    
    // Open input stream
    input.open(act_file);
    if (!input) {
        cerr << "Error opening file '" << act_file << "': " << strerror(errno) << "." << endl;
        exit(1);
    }

    // draw line from input
    unsigned int field_nr = 3; // initialize to 'valid' pattern
	while (getline(input, line)) {
        line_nr++;
		istringstream iss(line); // line as iss

        if (field_nr==3) field_nr=0;
        else {
            cerr << "Error! Line no. " << line_nr-1 << " in file '" << act_file << "' is not in correct input format. Please refer to README." << endl;
            exit(1);
        }

		while(getline(iss, tmp_field, '\t')) {	// split at tabs
			remove_dos_cr(&tmp_field);

			if (field_nr == 0) {		    // ID
                tid = (Tid) atoi(tmp_field.c_str());
                if (tid == 0) { 
                    cerr << "Error! Invalid ID: '" << tmp_field << "' in file '" << act_file << "', line " << line_nr << "." << endl; 
                    exit(1); 
                }
                field_nr=1;
			}

			else if (field_nr == 1) {	    // ACTIVITY NAME
				act_name = tmp_field;
                field_nr=2;
			}

			else if (field_nr == 2) {	    // ACTIVITY VALUES
                stringstream str;
                str  << tmp_field; int act_value; str >> act_value;
                if ((act_value != 0) && act_value != 1) { 
                    cerr << "Error! Invalid activity: '" << tmp_field << "' in file '" << act_file << "', line " << line_nr << "." << endl; 
                    exit(1);
                }
                fminer->AddActivity((bool) act_value, tid);
                field_nr=3;
			}

			else {
                cerr << "Error! Line no. " << line_nr << " in file '" << act_file << "' is not in correct input format. Please refer to README." << endl;
				exit(1);
			}

		}
	}


}


// main
int main(int argc, char *argv[], char *envp[]) {

    const char* program_name = argv[0];

    float def_chisq = 0.95;
    float chisq_sig = def_chisq;

    Frequency def_minfreq = 2;
    Frequency minfreq = def_minfreq;

    int def_type = 2;
    int type = def_type;
    
    int status=1;
    char* graph_file = NULL;
    char* act_file = NULL;

    bool do_output = true;
    bool refine_singles = false;
    bool aromatic = false;
    bool adjust_ub = true;
    bool do_pruning = true;
    bool do_backbone = true;
    bool line_nrs = false;
    bool bbrc_sep = false;
    bool most_specific_trees_only = false;


    // FILE ARGUMENT READ
    if (argc>2) {
       if (argv[argc-2][0]!='-') {
           graph_file = argv[argc-2]; status=0;
           if (argv[argc-1][0]=='-') {
               status=1;
           }
           else {
               act_file = argv[argc-1]; //chisq.active=1;
           }
       }
       else {
           if (argv[argc-1][0]=='-') {
               status=1;
           }
           else {
              graph_file = argv[argc-1]; //chisq.active=0;
              status = 0;
           }
       }
    }
    else if (argc==2){
       if (argv[argc-1][0]=='-') {
           status=1;
       }
       else {
           graph_file = argv[argc-1]; //chisq.active=0;
           status = 0;
       }
    }
    else status=1;

    
    // OPTIONS ARGUMENT READ
    char c;
    const char* const short_options = "f:l:p:saubmdonrh";
    const struct option long_options[] = {
        {"minfreq",                1, NULL, 'f'},
        {"level",                  1, NULL, 'l'},
        {"p-value",                1, NULL, 'p'},
        {"refine-singles",         0, NULL, 's'},
        {"no-aromaticity",         0, NULL, 'a'},
        {"no-upper-bound-pruning", 0, NULL, 'u'},
        {"no-bbr-classes",         0, NULL, 'b'},
        {"max-trees",              0, NULL, 'm'},
        {"no-dynamic-ub",          0, NULL, 'd'},
        {"no-output",              0, NULL, 'o'},
        {"line-nrs",               0, NULL, 'n'},
        {"bbrc-sep",               0, NULL, 'r'},
        {"help",                   0, NULL, 'h'},
        {NULL,                     0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
        switch(c) {
        case 'f':
            minfreq = atoi(optarg);
            break;
        case 'l':
            type = atoi (optarg);
            break;
        case 's':
            refine_singles = true;
            break;
        case 'a':
            aromatic = true;
            if (!graph_file) status = 1;
            break;
        case 'p':
            chisq_sig = atof (optarg);
            if (!act_file) status = 1;
            break;
        case 'u':
            do_pruning = false;
            if (!act_file) status = 1;
            break;
        case 'b':
            do_backbone = false;
            if (!act_file) status = 1;
            break;
        case 'm':
            most_specific_trees_only = true;
        case 'd':
            adjust_ub = false;
            if (!act_file) status = 1;
            break;
        case 'o':
            do_output = false;
            break;
        case 'n':
            line_nrs = true;
            break;
        case 'r':
            bbrc_sep = true;
            break;
        case 'h':
            status=2;
            break;
        case '?':
            status=2;
            break;
        default: 
            abort();
        }
    }


    // INTEGRITY CONSTRAINTS AND HELP OUTPUT
    if ((adjust_ub && !do_pruning) || (!do_backbone && adjust_ub) || (do_backbone && most_specific_trees_only)) status = 1;
             

    bool input_smi = false, input_gsp = false;
    string graph_file_str;
    if (graph_file) {
        graph_file_str = graph_file;
        string graph_file_suffix = graph_file_str.substr(graph_file_str.find_last_of("."));
        if (graph_file_suffix == ".smi") { input_smi=true; }
        else if (graph_file_suffix == ".gsp") { input_gsp=true; }
        else { cerr << "Suffix " << graph_file_suffix << " unknown!" << endl; status=2;}
    }

    if (status > 0) {
        cerr << "Usage: " << program_name << " [-f minfreq] [-l type] [-s] [-a] [-o] [-n] [-r] [-d [-b [-m] | -u]] [-p p_value] <graphs> <activities>" << endl;
        cerr << "       " << program_name << " [-f minfreq] [-l type] [-s] [-a] [-o] [-n] [-r] <graphs>" << endl;
        cerr << endl;
    }
    if (status==1) {
        cerr << "               use '-h' for additional information." << endl;
        return 1;
    }
    if (status > 1) {
        cerr << "  File formats:" << endl;
        cerr << "       <graphs> File must either have suffix .smi or .gsp, indicating SMILES or gSpan format." << endl;
        cerr << "       <activities> File must be in Activity format (suffix not relevant)." << endl;
        cerr << endl;
        cerr << "  General options:" << endl;
        cerr << "       -f  --minfreq _minfreq_      Set minimum frequency. Allowable values for _minfreq_: 1, 2, ... (default: " << def_minfreq<< ")." << endl;
        cerr << "       -l  --level _level_          Set fragment type. Allowable values for _type_: 1 (paths) and 2 (trees) (default: " << def_type << ")." << endl;
        cerr << "       -s  --refine-singles         Switch on refinement of fragments with frequency 1 (default: off)." << endl;
        cerr << "       -o  --no-output              Switch off output (default: on)." << endl;
        cerr << "       -n  --line-nrs               Switch on line numbers in output file (default: off)." << endl;
        cerr << "       -r  --bbrc-sep               Switch on BBRC separator in result vector (default: off)." << endl;
        cerr << endl;
        cerr << "  Upper bound pruning options:" << endl;
        cerr << "       -a  --aromaticity            Switch on aromatic ring perception when using smiles input format (default: off)." << endl;
        cerr << "       -d  --no-dynamic-ub          Switch off dynamic adjustment of upper bound for backbone mining (default: on)." << endl;
        cerr << "       -b  --no-bbr-classes         Switch off mining for backbone refinement classes (default: on)." << endl;
        cerr << "       -m  --max-trees              Switch on mining for maximal trees, aka the positive border (default: off)." << endl;
        cerr << "       -u  --no-upper-bound-pruning Switch off upper bound pruning (default: on)." << endl;
        cerr << "       -p  --p-value _p_value_      Set significance type. Allowable values for _p_value_: 0 <= _p_value_ <= 1.0 (default: " << def_chisq << ")." << endl;
        cerr << endl;
        cerr << "See README for additional information." << endl;
        cerr << endl;
        return 1;
    }  



    // DEFAULT SETTINGS FOR THIS HOST APP
   if (graph_file && act_file) {
        fminer = new Fminer(type, minfreq, chisq_sig, do_backbone);
        fminer->SetDynamicUpperBound(adjust_ub);
        fminer->SetPruning(do_pruning);
    }
    else if (graph_file) {
        fminer = new Fminer(type, minfreq);
        fminer->SetChisqActive(false);
    }
    else {
        exit(1);
    }
    fminer->SetAromatic(aromatic);
    fminer->SetRefineSingles(refine_singles);
    fminer->SetConsoleOut(true);
    fminer->SetDoOutput(do_output);
    fminer->SetLineNrs(line_nrs);
    fminer->SetBbrcSep(bbrc_sep);
    fminer->SetMostSpecTreesOnly(most_specific_trees_only);

    
    //////////
    // READ //
    //////////

    if (graph_file && act_file) {
        cerr << "Reading compounds..." << endl;
        if (input_smi) read_smi (graph_file);
        else if (input_gsp) read_gsp (graph_file);
        cerr << "Reading activities..." << endl;
        read_act (act_file);
    }
    
    else if (graph_file) {
        cerr << "Reading compounds..." << endl;
        if (input_smi) read_smi(graph_file);
        else if (input_gsp) read_gsp(graph_file);
    }

    //////////
    // MINE //
    //////////
    
    if (act_file) cerr << "Mining fragments... (bb: " << do_backbone << ", pr: " << do_pruning << ", adjub: " << adjust_ub << ", chisq sig: " << chisq_sig << ", min freq: " << minfreq << ", type: " << type << ")" << endl;
    else cerr << "Mining fragments... (min freq: " << minfreq << ", type: " << type << ")" << endl;

    clock_t t1 = clock ();
    for ( int j = 0; j < (int) fminer->GetNoRootNodes(); j++ ) {
        vector<string>* result = fminer->MineRoot(j);
        if (!fminer->GetConsoleOut()) { 
            each (*result) {
                cout << (*result)[i] << endl;
            }
        }
    }
    clock_t t2 = clock ();
//    if (!fminer->GetBackbone()) statistics->print();
    cerr << "Approximate total runtime: " << ( (float) t2 - t1 ) / CLOCKS_PER_SEC << "s" << endl;

    delete fminer;

}
