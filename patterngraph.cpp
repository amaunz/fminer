// patterngraph.cpp
// Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.
#include "patterngraph.h"
#include "graphstate.h"

void PatternGraph::init ( vector<CloseLegPtr> &closelegssource, int legindex ) {
  closetuples.push_back ( closelegssource[legindex]->tuple );
  frequency = closelegssource[legindex]->occurrences.frequency;
  
  // ADDED
  graphstate.closetuples = &closetuples;
  
  this->closelegssource = &closelegssource;
  this->legindex = legindex;
  
}


extern bool f;
extern string outl;

void PatternGraph::expand () {
  // ADDED
  int id = graphstate.is_normal ();
  if ( id == 0 ) {

    for ( int k = legindex + 1; k < (int) closelegssource->size (); k++ ) {
      if ( (*closelegssource)[k]->copy ) {
        CloseLegOccurrencesPtr closelegoccurrencesptr = 
          join ( (*closelegssource)[legindex]->occurrences, (*closelegssource)[k]->occurrences );
        if ( closelegoccurrencesptr ) {
          CloseLegPtr closelegptr = new CloseLeg;
          closelegptr->tuple = (*closelegssource)[k]->tuple;
          swap ( *closelegoccurrencesptr, closelegptr->occurrences );
          closelegs.push_back ( closelegptr );
        }
      }
    }

    // OUTPUT
    outl = graphstate.to_s(frequency);
    cout << outl;

    int addsize = statistics.patternsize + graphstate.edgessize - graphstate.nodes.size ();
    if ( addsize >= (int) statistics.frequenttreenumbers.size () ) {
      statistics.frequenttreenumbers.resize ( addsize + 1, 0 );
      statistics.frequentpathnumbers.resize ( addsize + 1, 0 );
      statistics.frequentgraphnumbers.resize ( addsize + 1, 0 );
    }
    statistics.frequentgraphnumbers[addsize]++;
    
    for ( int k = closelegs.size () - 1; k >= 0; k-- ) {

//      cerr << "Grow Graph" << endl;

      // GRAPHSTATE
      graphstate.insertEdge ( closelegs[k]->tuple.from, closelegs[k]->tuple.to, closelegs[k]->tuple.label );

      // TEST CHISQ HERE, BEFORE RECURSING!
      if (chisq.active) chisq.Calc(closelegs[k]->occurrences.elements);
//      cerr << "`-- " << setprecision(3) << chisq.p << " " << database.nodelabels[closelegs[k]->tuple.to].inputlabel << endl;

      // RECURSE
      PatternGraph patterngraph ( *this, k );
      patterngraph.expand ();
      graphstate.deleteEdge ( closelegs[k]->tuple.from, closelegs[k]->tuple.to );
    }
  }
  else 
    if ( id == 2 ) {
      (*closelegssource)[legindex]->copy = false; // should not be added to any later tree
    }
// cerr << "backtracking g" << endl; 
}

PatternGraph::~PatternGraph () {
  for ( int i = 0; i < (int) closelegs.size (); i++ )
    delete closelegs[i];
}
