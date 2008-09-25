// database.h
// Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.
#ifndef DATABASE_H
#define DATABASE_H
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include "legoccurrence.h"
#include "misc.h"

#include <openbabel/mol.h>
#include "openbabel/obconversion.h"


using namespace std;
using namespace OpenBabel;

typedef short InputEdgeLabel;
typedef string InputNodeLabel;
typedef short InputNodeId;
typedef unsigned int CombinedInputLabel;
//typedef OBAtom* OBAtomPtr;
//typedef OBBond* OBBondPtr;

#define combineInputLabels(label1,label2,label3) (label1 | ( ((unsigned int) label2 ) << 16 ) | ( ( (unsigned int) label3 ) << 24 ) )
// maximum 255 node labels for now.

#define NOINPUTEDGELABEL ((InputEdgeLabel) -1)
#define NOINPUTNODELABEL ((InputNodeLabel) -1)

template<class T>
class pvector {
public:
  T *array;
  int _size;
  pvector<T> ( T *array, int _size ): array ( array ), _size ( _size ) { }
  pvector<T> () { }
  inline int size () const { return _size; }
  void resize ( int s ) { _size = s; }
  void clear () { _size = 0; } // cannot remove allocation, as we are not managing that memory here 
  T &operator[] ( int i ) { return array[i]; }
};

struct DatabaseTreeEdge {
  EdgeLabel edgelabel;
  NodeId tonode;

  DatabaseTreeEdge ()  { }

  friend ostream &operator<< ( ostream &stream, DatabaseTreeEdge &databasetreeedge );
};

struct DatabaseTreeNode {
  NodeLabel nodelabel;
  bool incycle;
  pvector<DatabaseTreeEdge> edges;

  DatabaseTreeNode () { }

  friend ostream &operator<< ( ostream &stream, DatabaseTreeNode &databasetreenode );
};

/*
Nodes:			nodes =  [node1, ..., node n]
Edges			edges -> [e1 of n1,...,em of n1, ..., e1 of nn...,ek of nn]
*/

struct DatabaseTree {
  Tid tid, orig_tid;
  int line_nr;
  vector<DatabaseTreeNode> nodes;

  DatabaseTreeEdge *edges;
  int activity;

  DatabaseTree ( Tid tid , Tid orig_tid , int line_nr ): tid ( tid ), orig_tid (orig_tid ), line_nr (line_nr), activity ( -1 ) { }
  DatabaseTree () { }
  
  friend ostream &operator<< ( ostream &stream, DatabaseTree &databasetree );
};

typedef DatabaseTree *DatabaseTreePtr;

struct DatabaseNodeLabel {
  InputNodeLabel inputlabel;
  Frequency frequency;
  Tid lasttid;

  LegOccurrences occurrences;
  vector<EdgeLabel> frequentedgelabels;

  DatabaseNodeLabel (): frequency ( 1 ) { }
};

struct DatabaseEdgeLabel {
  InputEdgeLabel inputedgelabel;
  NodeLabel tonodelabel, fromnodelabel; 
  EdgeLabel edgelabel; // the (order) edge label to which this entry corresponds during the search
  Frequency frequency;
  Tid lasttid;

  DatabaseEdgeLabel (): frequency ( 1 ) { }
};

class Database {
  public:
    Database () { }
    vector<DatabaseTreePtr> trees;
    map<Tid, DatabaseTreePtr> trees_map;
    vector<DatabaseNodeLabel> nodelabels;
    vector<DatabaseEdgeLabel> edgelabels;
    map<InputNodeLabel,NodeLabel> nodelabelmap;
    map<CombinedInputLabel,EdgeLabel> edgelabelmap;
    vector<EdgeLabel> edgelabelsindexes; // given an edge label, returns the index of the element in edgelabels in which
    EdgeLabel frequentEdgeLabelSize () const { return edgelabelsindexes.size (); }
                                         // all information about this edge can be found. Used during the search,
					 // only frequent edge label, node label pairs are stored.

     // NOTE! In the input file, the nodes MUST be listed in pre-order.

    // reads a smiles string, and determines the frequency counts for node labels and (combined) edge labels
    void read_smi (char* smi_file);
    // reads binary activity information and stores it in the trees
    void read_act (char* smi_file);

     // after "read", determines the frequency of edges, using DatabaseNodeLabel's edgelasttid/edgelabelfrequency
    void edgecount ();

     // after "edgecount",
     // - removes infrequent data
     // - cleans up the datastructures used until now for counting frequencies
     // - changes the edge label order to optimise the search, fills the database with order numbers instead of
     //   the numbers assigned in the previous levels; fills edgelabelsindexes.
    void reorder ();

    void printTrees ();
    ~Database ();
  private:
    bool readTree (string smi, Tid tid , Tid orig_tid, int line_nr);
  
  	// Perform DFS through tree to identify cycles
    void determineCycledNodes ( DatabaseTreePtr tree, vector<int> &nodestack, vector<bool> &visited1, vector<bool> &visited2 );
};

extern Database database;

#endif
