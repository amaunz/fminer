// path.cpp
// Andreas Maunz, andreas@maunz.de, jul 2008
// Siegfried Nijssen, snijssen@liacs.nl, jan 2004.
#include <algorithm>
#include "patterntree.h"
#include "patterngraph.h"
#include "path.h"
#include "graphstate.h"
#include <iomanip>
#include <openbabel/data.h> 
#include "misc.h"

extern bool do_pruning;
extern bool do_backbone;
extern bool updated;
extern bool adjust_ub;
extern string outl;
extern int type;

// for every database node...
Path::Path ( NodeLabel startnodelabel ) {
  
    graphstate.insertStartNode ( startnodelabel );
    nodelabels.push_back ( startnodelabel );
    frontsymmetry = backsymmetry = totalsymmetry = 0;

    cerr << "Root: " << database.nodelabels[startnodelabel].inputlabel << endl; 

    DatabaseNodeLabel &databasenodelabel = database.nodelabels[startnodelabel];

    // ...gather frequent edge labels
    vector<EdgeLabel> frequentedgelabels;
    for ( unsigned int i = 0; i < databasenodelabel.frequentedgelabels.size (); i++ )
        frequentedgelabels.push_back ( database.edgelabels[databasenodelabel.frequentedgelabels[i]].edgelabel );
                                                                                                //  ^^^^^^^^^ is frequency rank!
    sort ( frequentedgelabels.begin (), frequentedgelabels.end () );                            // restores the rank order
    
    Tid lastself[frequentedgelabels.size ()];
    vector<EdgeLabel> edgelabelorder ( database.edgelabelsindexes.size () );
    EdgeLabel j = 0;

    // FOR ALL EDGES...
    for ( unsigned int i = 0; i < frequentedgelabels.size (); i++ ) {
        edgelabelorder[frequentedgelabels[i]] = j;                                              // For each rank, store it's position
        j++;
    
        // ...CREATE LEGS
        PathLegPtr leg = new PathLeg;
        legs.push_back ( leg );

        leg->tuple.depth = 0;                                           // TUPLE  DESCRIBES STRUCTURE...
        leg->tuple.connectingnode = 0;
        leg->tuple.edgelabel = frequentedgelabels[i]; 

        leg->occurrences.parent = &databasenodelabel.occurrences;       // ... OCCURRENCES DESCRIBES LOCATION IN TREE (1)
        leg->occurrences.number = 2;
        leg->occurrences.maxdegree = 0;
        leg->occurrences.selfjoin = 0;

        DatabaseEdgeLabel &databaseedgelabel = database.edgelabels[database.edgelabelsindexes[frequentedgelabels[i]]];
        leg->occurrences.frequency = databaseedgelabel.frequency;

        if ( databaseedgelabel.fromnodelabel == startnodelabel ) {
            leg->tuple.nodelabel = databaseedgelabel.tonodelabel;
        }

        else { 
            leg->tuple.nodelabel = databaseedgelabel.fromnodelabel;
        }

        lastself[i] = NOTID;

   }
    
    // ... OCCURRENCES DESCRIBES LOCATION IN TREE (2)
    for ( unsigned int i = 0; i < databasenodelabel.occurrences.elements.size (); i++ ) {
        DatabaseTree &tree = * (database.trees[databasenodelabel.occurrences.elements[i].tid]);
        DatabaseTreeNode &datanode = tree.nodes[databasenodelabel.occurrences.elements[i].tonodeid];
        for ( int j = 0; j < datanode.edges.size (); j++ ) {
            EdgeLabel edgelabel = edgelabelorder[datanode.edges[j].edgelabel];
            PathLeg &leg = * ( legs[edgelabel] );
            if ( !leg.occurrences.elements.empty () &&
                  leg.occurrences.elements.back ().occurrenceid == i &&
                  lastself[edgelabel] != tree.tid ) {
                leg.occurrences.selfjoin++;
                lastself[edgelabel] = tree.tid;
            }
            vector_push_back ( LegOccurrence, leg.occurrences.elements, legoccurrence );
            legoccurrence.tid = tree.tid;
            legoccurrence.occurrenceid = i;
            legoccurrence.tonodeid = datanode.edges[j].tonode;
            legoccurrence.fromnodeid = databasenodelabel.occurrences.elements[i].tonodeid;
        }
    }
  
}

Path::Path ( Path &parentpath, unsigned int legindex ) {
  PathLeg &leg = (*parentpath.legs[legindex]);
  int positionshift;
  
  // fill in normalisation information, it seems a lot of code, but in fact it's just a walk through the edge/nodelabels arrays.

  nodelabels.resize ( parentpath.nodelabels.size () + 1 );
  edgelabels.resize ( parentpath.edgelabels.size () + 1 );

  addCloseExtensions ( closelegs, parentpath.closelegs, leg.occurrences );

  if ( parentpath.nodelabels.size () == 1 ) {
    totalsymmetry = parentpath.nodelabels[0] - leg.tuple.nodelabel;
    frontsymmetry = backsymmetry = 0;
    nodelabels[1] = leg.tuple.nodelabel;
    edgelabels[0] = leg.tuple.edgelabel;
    nodelabels[0] = parentpath.nodelabels[0];
    positionshift = 0;
  }
  else if ( leg.tuple.depth == 0 ) {
    positionshift = 1;
    nodelabels[0] = leg.tuple.nodelabel;
    edgelabels[0] = leg.tuple.edgelabel;

    backsymmetry = parentpath.totalsymmetry;
    frontsymmetry = leg.tuple.nodelabel - parentpath.nodelabels[parentpath.nodelabels.size () - 2];
    totalsymmetry = leg.tuple.nodelabel - parentpath.nodelabels.back ();
    if ( !totalsymmetry )
      totalsymmetry = leg.tuple.edgelabel - parentpath.edgelabels.back ();

    unsigned int i = 0;
    // we can prepend only before strings of length 2
    if ( parentpath.nodelabels.size () > 2 ) {
      if ( !frontsymmetry )
        frontsymmetry = leg.tuple.edgelabel - parentpath.edgelabels[parentpath.nodelabels.size () - 3];

      while ( !frontsymmetry && i < parentpath.edgelabels.size () / 2 ) {
        nodelabels[i + 1] = parentpath.nodelabels[i];
        edgelabels[i + 1] = parentpath.edgelabels[i];

        frontsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 3];
        if ( !frontsymmetry && parentpath.nodelabels.size () > 3 )
          frontsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 4];

	    if ( !totalsymmetry ) {
	        totalsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 2];
	    if ( !totalsymmetry )
	        totalsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 3];
	    }

        i++;
      }
    }
    for ( ; !totalsymmetry && i < parentpath.edgelabels.size () / 2; i++ ) {
      nodelabels[i + 1] = parentpath.nodelabels[i];
      edgelabels[i + 1] = parentpath.edgelabels[i];

      totalsymmetry = parentpath.nodelabels[i] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 2];
      if ( !totalsymmetry && parentpath.nodelabels.size () > 3 )
        totalsymmetry = parentpath.edgelabels[i] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 3];
    }
    for ( ;i < parentpath.edgelabels.size (); i++ ) {
      nodelabels[i + 1] = parentpath.nodelabels[i];
      edgelabels[i + 1] = parentpath.edgelabels[i];
    }

    nodelabels[i + 1] = parentpath.nodelabels[i];



    // build OccurrenceLists
    extend ( leg.occurrences );
    for (unsigned int i = 0; i < candidatelegsoccurrences.size (); i++ ) {
      if ( candidatelegsoccurrences[i].frequency >= minfreq ) {
        PathLegPtr leg2 = new PathLeg;
        legs.push_back ( leg2 );
        leg2->tuple.edgelabel = i;
    	leg2->tuple.connectingnode = graphstate.lastNode ();
        DatabaseEdgeLabel &databaseedgelabel = database.edgelabels[database.edgelabelsindexes[i]];
        if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
          leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
        else
          leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
        leg2->tuple.depth = 0;
        store ( leg2->occurrences, candidatelegsoccurrences[i] ); // avoid copying
      }
    }



  }


  else  {
    positionshift = 0;

    frontsymmetry = parentpath.totalsymmetry;
    backsymmetry = parentpath.nodelabels[1] - leg.tuple.nodelabel;
    totalsymmetry = parentpath.nodelabels[0] - leg.tuple.nodelabel;
    if ( !totalsymmetry )
      totalsymmetry = parentpath.edgelabels[0] - leg.tuple.edgelabel;
    unsigned int i = 0;
    if ( parentpath.nodelabels.size () > 2 ) {
      if ( !backsymmetry )
        backsymmetry = parentpath.edgelabels[1] - leg.tuple.edgelabel;

      while ( !backsymmetry && i < parentpath.edgelabels.size () / 2 ) {
        nodelabels[i] = parentpath.nodelabels[i];
	edgelabels[i] = parentpath.edgelabels[i];

	backsymmetry = parentpath.nodelabels[i + 2] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
	if ( !backsymmetry && parentpath.nodelabels.size () > 3 )
	  backsymmetry = parentpath.edgelabels[i + 2] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];

	if ( !totalsymmetry ) {
	  totalsymmetry = parentpath.nodelabels[i + 1] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
	  if ( !totalsymmetry && parentpath.nodelabels.size () > 3 )
	    totalsymmetry = parentpath.edgelabels[i + 1] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];
	}
	i++;
      }
    }
    for ( ; !totalsymmetry && i < parentpath.edgelabels.size () / 2; i++ ) {
      nodelabels[i] = parentpath.nodelabels[i];
      edgelabels[i] = parentpath.edgelabels[i];
      totalsymmetry = parentpath.nodelabels[i + 1] - parentpath.nodelabels[parentpath.nodelabels.size () - i - 1];
      if ( !totalsymmetry && i < parentpath.edgelabels.size () - 1 )
        totalsymmetry = parentpath.edgelabels[i + 1] - parentpath.edgelabels[parentpath.nodelabels.size () - i - 2];
    }
    for ( ; i < parentpath.edgelabels.size (); i++ ) {
      nodelabels[i] = parentpath.nodelabels[i];
      edgelabels[i] = parentpath.edgelabels[i];
    }

    nodelabels[i] = parentpath.nodelabels[i];
    edgelabels[i] = leg.tuple.edgelabel;
    nodelabels[i+1] = leg.tuple.nodelabel;
  }

  unsigned int i = 0;
  LegOccurrencesPtr legoccurrencesptr;
  for ( ; i < legindex; i++ ) {
    PathLeg &leg2 = (*parentpath.legs[i]);

    if ( (legoccurrencesptr = join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) { // JOIN OCCURRENCES
      PathLegPtr leg3 = new PathLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( (legoccurrencesptr = join ( leg.occurrences )) ) {
    PathLegPtr leg3 = new PathLeg;
    legs.push_back ( leg3 );
    leg3->tuple.connectingnode = leg.tuple.connectingnode;
    leg3->tuple.edgelabel = leg.tuple.edgelabel;
    leg3->tuple.nodelabel = leg.tuple.nodelabel;
    leg3->tuple.depth = leg.tuple.depth + positionshift;
    store ( leg3->occurrences, *legoccurrencesptr );
  }

  for ( i++; i < parentpath.legs.size (); i++ ) {
    PathLeg &leg2 = (*parentpath.legs[i]);
    if ( (legoccurrencesptr = join ( leg.occurrences, leg2.tuple.connectingnode, leg2.occurrences )) ) {
      PathLegPtr leg3 = new PathLeg;
      legs.push_back ( leg3 );
      leg3->tuple.connectingnode = leg2.tuple.connectingnode;
      leg3->tuple.edgelabel = leg2.tuple.edgelabel;
      leg3->tuple.nodelabel = leg2.tuple.nodelabel;
      leg3->tuple.depth = leg2.tuple.depth + positionshift;
      store ( leg3->occurrences, *legoccurrencesptr );
    }
  }

  if ( positionshift ) {
    addCloseExtensions ( closelegs, leg.occurrences.number ); // stored separately
    return;
  }

  extend ( leg.occurrences );
  for ( unsigned int i = 0; i < candidatelegsoccurrences.size (); i++ ) {
    if ( candidatelegsoccurrences[i].frequency >= minfreq ) {
      PathLegPtr leg2 = new PathLeg;
      legs.push_back ( leg2 );
      leg2->tuple.edgelabel = i;
      leg2->tuple.connectingnode = graphstate.lastNode ();
      DatabaseEdgeLabel &databaseedgelabel = database.edgelabels[database.edgelabelsindexes[i]];
      if ( databaseedgelabel.fromnodelabel == leg.tuple.nodelabel )
        leg2->tuple.nodelabel = databaseedgelabel.tonodelabel;
      else
        leg2->tuple.nodelabel = databaseedgelabel.fromnodelabel;
      leg2->tuple.depth = leg.tuple.depth + 1;
      store ( leg2->occurrences, candidatelegsoccurrences[i] ); // avoid copying
    }
  }

  addCloseExtensions ( closelegs, leg.occurrences.number );
}

Path::~Path () {
  for ( unsigned int i = 0; i < legs.size (); i++ )
    delete legs[i];
  for ( unsigned int i = 0; i < closelegs.size (); i++ )
    delete closelegs[i];
}

// ADDED
bool Path::is_normal ( EdgeLabel edgelabel ) {
  // symplistic quadratic algorithm
  int nodelabelssize = nodelabels.size (), step, add, start;
  
  edgelabels.push_back ( edgelabel );
    
  // if we would program it better, we would use the 'totalsymmetry' variable here;
  // however, to be quick & easy, we used a different coding here...
  int t = nodelabelssize - 1, r = 0;
  int symmetry;
  do {
    symmetry = nodelabels[t] - nodelabels[r];
    int nt = ( t + nodelabelssize - 1 ) % nodelabelssize;
    if ( !symmetry ) 
      symmetry = edgelabels[nt] - edgelabels[r];
    r = ( r + 1 ) % nodelabelssize;
    t = nt;
  }
  while ( symmetry == 0 && t != nodelabelssize - 1 );
  
  if ( symmetry < 0 ) {
    step = -1 + nodelabelssize ;
    add = -1 + nodelabelssize ;
    start = nodelabelssize - 1;
  }
  else {
    step = 1 + nodelabelssize;
    add = nodelabelssize ;
    start = 0;
  }
  for ( int i = 0; i < nodelabelssize; i++ ) {
    // starting positions for the new path
    int k = start, l = i, p;
    do {
      if ( nodelabels[l] < nodelabels[k] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( nodelabels[l] > nodelabels[k] )
        break;
      p = ( k + add ) % nodelabelssize;
      l = ( l + nodelabelssize - 1 ) % nodelabelssize;
      if ( edgelabels[l] < edgelabels[p] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( edgelabels[l] > edgelabels[p] ) 
        break;
      k = ( k + step ) % nodelabelssize;
    }
    while ( k != start );
    
    k = start, l = i;
    do {
      if ( nodelabels[l] < nodelabels[k] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( nodelabels[l] > nodelabels[k] ) 
        break;
      p = ( k + add ) % nodelabelssize;
      if ( edgelabels[l] < edgelabels[p] ) {
        edgelabels.pop_back ();
        return false;
      }
      if ( edgelabels[l] > edgelabels[p] ) 
        break;
      l = ( l + 1 ) % nodelabelssize;
      k = ( k + step ) % nodelabelssize;
    }
    while ( k != start );
    
  }
  edgelabels.pop_back ();
  return true;
}









void Path::expand2 (pair<float,string> max) {

  vector<unsigned int> forwpathlegs; forwpathlegs.clear();
  vector<unsigned int> backwpathlegs; backwpathlegs.clear();
  vector<unsigned int> pathlegs; pathlegs.clear();

  for ( unsigned int i = 0; i < legs.size (); i++ ) {

    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.depth == nodelabels.size () - 1 ) {

      if ( tuple.nodelabel > nodelabels[0] ||
           ( tuple.nodelabel == nodelabels[0] &&
             ( tuple.edgelabel > edgelabels[0] ||
               ( tuple.edgelabel == edgelabels[0] && backsymmetry <= 0 )
             )
           ) ) {
        forwpathlegs.push_back(i);
        pathlegs.push_back(i);

      }
    }
  }

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.depth != nodelabels.size () - 1 ) {


      if ( legs[i]->tuple.depth == 0 ) {
        if ( totalsymmetry &&
             ( tuple.nodelabel > nodelabels.back () ||
             ( tuple.nodelabel == nodelabels.back () &&
               ( tuple.edgelabel > edgelabels.back () ||
                 ( tuple.edgelabel == edgelabels.back () && frontsymmetry >= 0 )
               )
             ) ) ) {
            backwpathlegs.push_back(i);
            pathlegs.push_back(i);
        }
      }
    }
  }
 

  // FREE STRUCTURES: we have reached a leaf
  if (do_backbone && (pathlegs.size()==0)) { 
    if (updated) { 
        cout << max.second;
        updated = false;
    }
  }


  statistics.patternsize++;
  if ( (unsigned) statistics.patternsize > statistics.frequenttreenumbers.size () ) {
    statistics.frequenttreenumbers.push_back ( 0 );
    statistics.frequentpathnumbers.push_back ( 0 );
    statistics.frequentgraphnumbers.push_back ( 0 );
  }
  ++statistics.frequentpathnumbers[statistics.patternsize-1];
  
  if ( statistics.patternsize == maxsize ) {
    statistics.patternsize--;
    return;
  }


  
  // PHASE 3: INSERT CYLE CLOSING EDGES IN CASE OF CYCLE-CLOSING REFINEMENT
  if ( closelegs.size () && type > 2 ) {
    
    NodeId from = graphstate.nodes.size ();
    NodeId to = 0;
    while ( graphstate.nodes[to].edges.size () == 2 )
      to++;
    to++;
      
    for ( unsigned int i = 0; i < closelegs.size (); i++ ) {
      if ( closelegs[i]->tuple.from == from &&
           closelegs[i]->tuple.to == to &&
           is_normal ( closelegs[i]->tuple.label ) ) {

           // Processing Closelegs

           // GRAPHSTATE
           graphstate.insertEdge ( closelegs[i]->tuple.from, closelegs[i]->tuple.to, closelegs[i]->tuple.label );

           // DO NOT RECURSE

           // OUTPUT
           if (chisq.active) chisq.Calc(closelegs[i]->occurrences.elements);
           outl = graphstate.to_s(closelegs[i]->occurrences.frequency);
           cout << outl;

           int addsize = statistics.patternsize + graphstate.edgessize - graphstate.nodes.size ();
           if ( (unsigned) addsize >= statistics.frequenttreenumbers.size () ) {
              statistics.frequenttreenumbers.resize ( addsize + 1, 0 );
              statistics.frequentpathnumbers.resize ( addsize + 1, 0 );
              statistics.frequentgraphnumbers.resize ( addsize + 1, 0 );
           }
           statistics.frequentgraphnumbers[addsize]++;
           graphstate.deleteEdge ( closelegs[i]->tuple.from, closelegs[i]->tuple.to );
        
           // DON'T RECURSE GRAPH GROWING!
           // circle graphs can only grow from paths, all other graphs
           // can grow from spanning trees, which we prefer for now
      }
    }
  }


  // Grow Path forw
  for (unsigned int j=0; j<forwpathlegs.size() ; j++ ) {
    unsigned int index = forwpathlegs[j];


    // Calculate chisq
    if (chisq.active) chisq.Calc(legs[index]->occurrences.elements);
          
    // GRAPHSTATE AND OUTPUT
    graphstate.insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );
    outl = graphstate.to_s(legs[index]->occurrences.frequency);

    // immediate output
    if (!do_backbone) cout << outl;

    // RECURSE
    float cmax = maxi ( maxi ( chisq.sig, max.first ), chisq.p );
    if ( (
             !do_pruning || 
             (
               (  !adjust_ub && (chisq.u >= chisq.sig) ) || 
               (   adjust_ub && (chisq.u >= cmax) )
             )
         ) &&
         (
            legs[index]->occurrences.frequency>1
         )
      ){   // UB-PRUNING

      Path path ( *this, index );
      if (max.first<chisq.p) { updated = true; path.expand2 ( pair<float, string>(chisq.p, outl)); }
      else path.expand2 (max);
    }
    else {
        if (do_backbone && updated) {  // FREE STRUCTURES: search was pruned
            cout << max.second;
            updated=false;
        }
    }

    graphstate.deleteNode ();
  }



  // Grow Path backw
  for (unsigned int j=0; j<backwpathlegs.size() ; j++ ) {
    unsigned int index = backwpathlegs[j];
    

    // Calculate chisq
    if (chisq.active) chisq.Calc(legs[index]->occurrences.elements);

    // GRAPHSTATE AND OUTPUT
    graphstate.insertNode ( legs[index]->tuple.connectingnode, legs[index]->tuple.edgelabel, legs[index]->occurrences.maxdegree );
    outl = graphstate.to_s(legs[index]->occurrences.frequency);

    // immediate output
    if (!do_backbone) cout << outl;

    // RECURSE
    float cmax = maxi ( maxi ( chisq.sig, max.first ), chisq.p );
    if ( ( 
             !do_pruning || 
             (
               (  !adjust_ub && (chisq.u >= chisq.sig) ) || 
               (   adjust_ub && (chisq.u >= cmax) )
             )
         ) &&
         (
            legs[index]->occurrences.frequency>1
         )
     ){   // UB-PRUNING

      Path path ( *this, index );
      if (max.first<chisq.p) { updated = true; path.expand2 ( pair<float, string>(chisq.p, outl)); }
      else path.expand2 (max);
    }
    else {
        if (do_backbone && updated) { // FREE STRUCTURES: search was pruned
            cout << max.second;
            updated=false;
        }
    }

    graphstate.deleteNode ();

  }



  

  bool uptmp = updated;

  
//  cerr << "BEGIN" << endl;
//  cin.get();

  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.depth != nodelabels.size () - 1 ) {


      // PHASE 2: GROW TREE
      if ( legs[i]->tuple.depth != 0 ) {
        if ( ( totalsymmetry || legs[i]->tuple.depth <= edgelabels.size () / 2 ) &&
 	      ( legs[i]->tuple.depth != 1 || legs[i]->tuple.edgelabel >= edgelabels[0] ) &&
	      ( legs[i]->tuple.depth != nodelabels.size () - 2 || legs[i]->tuple.edgelabel >= edgelabels.back () ) &&
	        type > 1 ) {
          // Calculate chisq
          if (chisq.active) chisq.Calc(legs[i]->occurrences.elements);
          // GRAPHSTATE
          graphstate.insertNode ( legs[i]->tuple.connectingnode, legs[i]->tuple.edgelabel, legs[i]->occurrences.maxdegree );
          outl = graphstate.to_s(legs[i]->occurrences.frequency);

          if (!do_backbone) cout << outl;

          // RECURSE
          float cmax = maxi ( maxi ( chisq.sig, max.first ), chisq.p );

          if ( ( !do_pruning || 
               (
                 (  !adjust_ub && (chisq.u >= chisq.sig) ) || 
                 (   adjust_ub && (chisq.u >= cmax) )
               )
             ) &&
             (
                legs[i]->occurrences.frequency>1
             )
    
          ){   // UB-PRUNING

            PatternTree tree ( *this, i );
            if (max.first<chisq.p) { updated = true; tree.expand ( pair<float, string>(chisq.p, outl) ); }
            else tree.expand (max);
          }

          else {
            if (do_backbone && updated) { 
              cout << max.second;
              updated=false;
            }
          }

	      graphstate.deleteNode ();
        }
      }
    }
  }

//  cerr << "END" << endl;
//  cin.get();



  updated=uptmp;
    
  statistics.patternsize--;


//  cerr << "backtracking p" << endl;
}








void Path::expand () {


  for ( unsigned int i = 0; i < legs.size (); i++ ) {
    PathTuple &tuple = legs[i]->tuple;
    if ( tuple.nodelabel >= nodelabels[0] ) {
        
      if (chisq.active) chisq.Calc(legs[i]->occurrences.elements);

      // GRAPHSTATE AND OUTPUT
      graphstate.insertNode ( tuple.connectingnode, tuple.edgelabel, legs[i]->occurrences.maxdegree );
      outl = graphstate.to_s(legs[i]->occurrences.frequency);
      if (!do_backbone) cout << outl;

      // RECURSE
      Path path (*this, i);
      updated = true;
      path.expand2 (pair<float, string>(chisq.p, outl));
      graphstate.deleteNode ();

    }
  }
  graphstate.deleteStartNode ();


//  cerr << "backtracking p" << endl;
}






ostream &operator<< ( ostream &stream, Path &path ) {
  stream << /* database.nodelabels[ */ (int) path.nodelabels[0] /* ].inputlabel; */ << " ";
  for ( unsigned int i = 0; i < path.edgelabels.size (); i++ ) {
    //stream << (char) ( path.edgelabels[i] + 'A' ) << path.nodelabels[i+1];
    stream << /*database.edgelabels[database.edgelabelsindexes[*/ (int) path.edgelabels[i] /*]].inputedgelabel */ << " " <<  /* database.nodelabels[ */ (int) path.nodelabels[i+1] /* ].inputlabel */ << " ";
  }
  return stream;
}
