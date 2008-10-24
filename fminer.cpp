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

#include <openbabel/mol.h>
#include "openbabel/obconversion.h"

#include <getopt.h>
#include <time.h>
#include <iostream>
#include <string.h>
#include <sstream>

#include "database.h"
#include "misc.h"
#include "legoccurrence.h"
#include "graphstate.h"
#include "constraints.h"

using namespace std;
using namespace OpenBabel;

// global constraints

float def_chisq = 0.95;
ChisqConstraint chisq(def_chisq);

Frequency def_minfreq = 2;
Frequency minfreq = def_minfreq;

int def_type = 2;
int type = def_type;

bool adjust_ub = true;
bool do_pruning = true;
bool do_backbone = true;
bool updated = true;

Database database; 
Statistics statistics; 
int maxsize = ( 1 << ( sizeof(NodeId)*8 ) ) - 1; // safe default for the largest allowed pattern
string outl = "";

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

void read_smi (char* smi_file) {
    Tid tree_nr = 0; 
    Tid tree_id = 0;
    int line_nr = 1;

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
                if (tree_id == 0) { cerr << "Error! Invalid ID: '" << tmp_field << "' in file " << smi_file << ", line " << line_nr << "." << endl; exit(1); }
            }
            else if (field_nr == 1) { // SMILES
                if (database.readTree (tmp_field , tree_nr, tree_id, line_nr)) {
                    tree_nr++;
                }
                else {
                    cerr << "on line " << line_nr << ", id " << tree_id << "." << endl;
                }
            }
            field_nr++;
        }
        line_nr++;
    }

    cerr << tree_nr << " compounds" << endl;
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
					if (database.trees_map[tid] == NULL) {	// ignore compounds without structures
						no_id = tmp_field;
						cerr << "No structure for ID " << tid << ". Ignoring entry!" << endl;
						break;
					}
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
                if ((database.trees_map[tid]->activity = (bool) act_value)) { chisq.na++; }
                else { chisq.ni++; }
			}
			else {
				cerr << "Error! More than 3 columns at line " << line_nr << "." << endl;
				exit(1);
			}
			field_nr++;
		}
		line_nr++;
	}

    chisq.n = chisq.na + chisq.ni;
}


// main
int main(int argc, char *argv[], char *envp) {

    
    int status=1;
    char* smi_file = NULL;
    char* act_file = NULL;

    if (argc>2) {
       if (argv[argc-2][0]!='-') {
           smi_file = argv[argc-2]; status=0;
           if (argv[argc-1][0]=='-') {
               status=1;
           }
           else {
               act_file = argv[argc-1]; chisq.active=1;
           }
       }
       else {
            status=1;
       }
    }
    else status=1;

    char c;

    while ((c = getopt(argc, argv, "p:l:f:cxjh")) != -1) {
        switch(c) {
        case 'p':
            chisq.sig = atof (optarg);
            if (chisq.sig < 0.0 || chisq.sig > 1.0) { cerr << "Error! Invalid value '" << chisq.sig << "' for parameter p." << endl; status = 2; }
            if (!act_file) status = 2;
            break;
        case 'l':
            type = atoi (optarg);
            if ((type != 1) && (type != 2)) { cerr << "Error! Invalid value '" << type << "' for parameter l." << endl; status = 2; }
            break;
        case 'f':
            minfreq = atoi(optarg);
            if (minfreq < 1) { cerr << "Error! Invalid value '" << minfreq << "' for parameter f." << endl; status = 2; }
            break;
        case 'c':
            do_backbone = false;
            if (!chisq.active) status = 2;
            break;
        case 'x':
            do_pruning = false;
            if (!chisq.active) status = 2;
            break;
        case 'j':
            adjust_ub = false;
            if (!do_backbone) status = 2;
            if (!do_pruning) status = 2;
            break;
       case 'h':
            status=2;
            break;

        default: abort();
        }
    }

    if (status > 0) {
        cerr << "fminer: usage: fminer [-f minfreq] [-l type] [-p p_value] [ [-x] [-c] | [-j] ] [-h] smiles_file activities_file" << endl;
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
        cerr << endl;
        cerr << "See README for additional information." << endl;
        cerr << endl;
        return 1;
    }  



    //////////
    // READ //
    //////////

    cerr << "Reading compounds..." << endl;
    read_smi (smi_file);
    
    if (chisq.active) {
        cerr << "Reading activities..." << endl;
        read_act (act_file);
        each (database.trees) {
            if (database.trees[i]->activity == -1) {
                cerr << "Error! ID " << database.trees[i]->orig_tid << " is missing activity information." << endl;
                exit(1);
            }
        }
    }
    else do_pruning = false; // ensure every recursion is done

    database.edgecount ();
    database.reorder ();



    //////////
    // MINE //
    //////////
    

    if (!do_pruning || !do_backbone) adjust_ub = false; // ensure adjust_ub only for pruning and backbone
    if (!act_file) do_backbone = false;

    cerr << "Mining fragments... (bb: " << do_backbone << ", pr: " << do_pruning << ", adjub: " << adjust_ub << ", chisq sig: " << chisq.sig << ", min freq: " << minfreq << ")" << endl;
    
    initLegStatics ();
    graphstate.init ();

//    if (chisq.active) cout << "# - [ smiles,    p_chisq,    occ_list_active,    occ_list_inactive ]\n";
//    else              cout << "# - [ smiles,    frequency ]\n";
 
    clock_t t1 = clock ();
    for ( int i = 0; i < (int) database.nodelabels.size (); i++ ) {
        if ( database.nodelabels[i].frequency >= minfreq && database.nodelabels[i].frequentedgelabels.size () ) {
            Path path(i);
            path.expand ();
        }
    }
    clock_t t2 = clock ();

    cerr << "Approximate total runtime: " << ( (float) t2 - t1 ) / CLOCKS_PER_SEC << "s" << endl;

}
