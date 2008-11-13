/* Copyright (C) 2008  Andreas Maunz <andreas@maunz.de> 
 *
 *    
 *        This program is free software; you can redistribute it and/or modify
 *        it under the terms of the GNU General Public License as published by
 *        the Free Software Foundation; either version 2 of the License, or
 *        (at your option) any later version.
 *
 *        This program is distributed in the hope that it will be useful,
 *        but WITHOUT ANY WARRANTY; without even the implied warranty of
 *        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *        GNU General Public License for more details.
 *
 *        You should have received a copy of the GNU General Public License
 *        along with this program; if not, write to the Free Software
 *        Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <getopt.h>
#include <time.h>
#include <iostream>
#include <string.h>

#include "fminer.h"

using namespace std;
using namespace OpenBabel;

extern Statistics* statistics;

Fminer* fm;

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

void read_gsp (char* gsp_file) {
    FILE *input = fopen (gsp_file, "r");
    fm->ReadGsp(input);
}

void read_smi (char* smi_file) {
    Tid tree_id = 0;

    ifstream input;
    string line;
    string tmp_field;
    input.open(smi_file);

    if (!input) {
        cerr << "Cannot open " << smi_file << endl;
        exit(1);
    }

    while (getline(input, line)) {
        istringstream iss(line);
        string id = "";
        int field_nr = 0;
        while(getline(iss, tmp_field, '\t')) {  // split at tabs
            if (field_nr == 0) { // ID
                tree_id = (unsigned int) atoi (tmp_field.c_str());
            }
            else if (field_nr == 1) { // SMILES
                fm->AddCompound (tmp_field , tree_id);

            }
            field_nr++;

        }
    }

    cerr << fm->GetNoCompounds() << " compounds" << endl;
    input.close();

}

void read_act (char* act_file) {
    ifstream input;
    string line;
    string tmp_field;
    string no_id;
    string act_name;
    Tid tid=0;
    unsigned int line_nr = 0;
    
    input.open(act_file);
    if (!input) {
        cerr << "Cannot open " << act_file << endl;
        exit(1);
    }

	while (getline(input, line)) {
		istringstream iss(line); 
		unsigned int field_nr = 0;
		while(getline(iss, tmp_field, '\t')) {	// split at tabs
			remove_dos_cr(&tmp_field);
			if (field_nr == 0) {		// ID
				if (tmp_field == no_id)		// ignore compounds without structures
					break;
				else {
                    tid = (Tid) atoi(tmp_field.c_str());
                    if (tid == 0) { cerr << "Error! Invalid ID: '" << tmp_field << "' in file " << act_file << ", line " << line_nr+1 << "." << endl; exit(1); }
				}
			}
			else if (field_nr == 1) {	// ACTIVITY NAME
				act_name = tmp_field;
			}
			else if (field_nr == 2) {	// ACTIVITY VALUES
                stringstream str;
                str  << tmp_field;
                int act_value;
                str >> act_value;
                if ((act_value != 0) && act_value != 1) { cerr << "Error! Invalid activity: '" << tmp_field << "' in file " << act_file << ", line " << line_nr+1 << "." << endl; exit(1); }
                fm->AddActivity((bool) act_value, tid);
			}
			else {
				cerr << "Error! More than 3 columns at line " << line_nr << "." << endl;
				exit(1);
			}
			field_nr++;
		}
		line_nr++;
	}


}


// main
int main(int argc, char *argv[], char *envp) {


    float def_chisq = 0.95;
    float chisq_sig = def_chisq;

    Frequency def_minfreq = 2;
    Frequency minfreq = def_minfreq;

    int def_type = 2;
    int type = def_type;
    
    int status=1;
    char* smi_file = NULL;
    char* act_file = NULL;
    char* gsp_file = NULL;

    bool refine_singles = false;
    bool aromatic = true;
    bool adjust_ub = true;
    bool do_pruning = true;
    bool do_backbone = true;



    // FILE ARGUMENT READ
    if (argc>2) {
       if (argv[argc-2][0]!='-') {
           smi_file = argv[argc-2]; status=0;
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
              gsp_file = argv[argc-1]; //chisq.active=0;
              status = 0;
           }
       }
    }
    else if (argc==2){
       if (argv[argc-1][0]=='-') {
           status=1;
       }
       else {
           gsp_file = argv[argc-1]; //chisq.active=0;
           status = 0;
       }
    }
    else status=1;


    
    // OPTIONS ARGUMENT READ
    char c;
    while ((c = getopt(argc, argv, "p:l:f:cxjhas")) != -1) {
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
            aromatic = false;
            if (!smi_file) status = 1;
            break;
        case 'p':
            chisq_sig = atof (optarg);
            if (gsp_file) status = 1;
            break;
        case 'x':
            do_pruning = false;
            if (gsp_file) status = 1;
            break;
        case 'c':
            do_backbone = false;
            if (gsp_file) status = 1;
            break;
       case 'j':
            adjust_ub = false;
            if (gsp_file) status = 1;
            break;
       case 'h':
            status=2;
            break;
       default: abort();
        }
    }


    // INTEGRITY CONSTRAINTS AND HELP OUTPUT
    if ((!do_pruning && !do_backbone) || (!do_backbone && !adjust_ub) || (!adjust_ub) && (!do_pruning)) status = 1;


    if (status > 0) {
        cerr << "fminer usage: fminer [-f minfreq] [-l type] [-s] [-a] [-p p_value] [ -x | -c | -j ] smiles_file activities_file" << endl;
        cerr << "              fminer [-f minfreq] [-l type] [-s] gspan_file" << endl;
    }
    if (status==1) {
        cerr << "               use '-h' for additional information." << endl;
        return 1;
    }
    if (status > 1) {
        cerr << "       -f Set minimum frequency. Allowable values for _minfreq_: 1, 2, ... Default is " << def_minfreq<< "." << endl;
        cerr << "       -l Set fragment type. Allowable values for _type_: 1 (paths) and 2 (trees). Default is " << def_type << "." << endl;
        cerr << "       -p Set significance type. Allowable values for _p_value_: 0 <= p_value <= 1.0. Default is " << def_chisq << "." << endl;
        cerr << "       -x Switch off upper bound pruning (less efficient, use only for performance evaluation)." << endl;
        cerr << "       -c Switch off backbone mining." << endl;
        cerr << "       -j Switch off dynamic adjustment of upper bound for backbone mining (less efficient, use only for performance evaluation). Implied by both -x or -c." << endl;
        cerr << "       -s Switch on refinement of fragments with frequency 1 (default: off)." << endl;
        cerr << "       -a Switch off aromatic ring perception when using smiles input format (default: on)." << endl;
        cerr << endl;
        cerr << "See README for additional information." << endl;
        cerr << endl;
        return 1;
    }  



    // DEFAULT SETTINGS FOR THIS HOST APP
   if (smi_file && act_file) {
        fm = new Fminer(type, minfreq, chisq_sig, do_backbone);
        fm->SetDynamicUpperBound(adjust_ub);
        fm->SetPruning(do_pruning);
    }
    else if (gsp_file) {
        fm = new Fminer(type, minfreq);
        fm->SetChisqActive(false);
        fm->SetDynamicUpperBound(false);
        fm->SetPruning(false);
        fm->SetBackbone(false);
 
    }
    else {
        exit(1);
    }
    fm->SetAromatic(aromatic);
    fm->SetRefineSingles(refine_singles);
    fm->SetConsoleOut(true);
 

    
    //////////
    // READ //
    //////////

    if (smi_file && act_file) {
        cerr << "Reading compounds..." << endl;
        read_smi (smi_file);
        cerr << "Reading activities..." << endl;
        read_act (act_file);
    }
    
    else if (gsp_file) {
        read_gsp(gsp_file);
    }

    //////////
    // MINE //
    //////////
    
    if (!gsp_file) cerr << "Mining fragments... (bb: " << do_backbone << ", pr: " << do_pruning << ", adjub: " << adjust_ub << ", chisq sig: " << chisq_sig << ", min freq: " << minfreq << ", type: " << type << ")" << endl;
    else cerr << "Mining fragments... (min freq: " << minfreq << ", type" << type << ")" << endl;

    clock_t t1 = clock ();
    for ( int j = 0; j < (int) fm->GetNoRootNodes(); j++ ) {
        vector<string>* result = fm->MineRoot(j);
        if (!fm->GetConsoleOut()) { 
            each (*result) {
                cout << (*result)[i] << endl;
            }
        }
    }
    clock_t t2 = clock ();
    if (gsp_file) statistics->print();
    cerr << "Approximate total runtime: " << ( (float) t2 - t1 ) / CLOCKS_PER_SEC << "s" << endl;

    delete fm;

}
