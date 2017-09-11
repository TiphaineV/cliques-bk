#include <stack>
#include <string>
#include <sstream>
#include<set>
#include<vector>
#include<map>
#include<iterator>
#include<iostream>
#include<stdio.h>
#include<assert.h>
#include<stdlib.h>

#define MAX_LINE_LENGTH 1000



//clock_t is_clique_dicho_search = 0;
//clock_t map_find = 0;
//clock_t loopX_is_clique = 0;
//clock_t run_in_is_clique = 0;

typedef int node;


class interval{
    public:
        float b;
        float e;
};

typedef std::map<node, std::vector<interval> > linktype;

class nodeint{
    public:
        node v;
        float b;
        float e;
};

class Clique{
    public :
        float b;
        float e;
        std::pair< std::set<node>, std::set<node> > X;

        bool operator <( const Clique &rhs ) const
        {
            if( b != rhs.b )
                return( b < rhs.b );
            if( e != rhs.e )
                return( e < rhs.e);
            int sX_top = X.first.size();
            int sX_bot = X.second.size();
            int srhsX_top = rhs.X.first.size();
            int srhsX_bot = rhs.X.second.size();
            if( sX_top != srhsX_top )
                return( sX_top < srhsX_top );
            if( sX_bot != srhsX_bot )
                return( sX_top < srhsX_bot );
            //return( X < rhs.X );
            // Comment comparer la taille de 2 cliques biparties ? |top|+|bot| ?
            return( (sX_top + sX_bot) < (srhsX_top + srhsX_bot) );
        }

        //Constructor?
        Clique(float bc, float ec, std::pair<std::set<node>, std::set<node> > Xc);

        //Destructor?
        ~Clique();

        void print_clique(FILE *F) const ;
};

Clique::Clique(float bc, float ec, std::pair< std::set<node>, std::set<node> > Xc){
    //fprintf(stdout, "%s\n", "Constructing clique");
    //fprintf(stdout, "%p\n", &X);
    b = bc;
    e = ec;
    X.first = std::set<node>(Xc.first);
    X.second = std::set<node>(Xc.second);
}

Clique::~Clique(){
    //fprintf(stdout, "%s\n", "Calling ~Clique");
    //X.~set();
}


void Clique::print_clique(FILE *F) const {
    fprintf(F, "%f %f ", b, e);
    //debug
    std::pair<int, int> Xsize = std::pair<int,int>(0,0);

    // print nodes
    // std::cout << ", nodes : ";
    for( std::set<node>::iterator itX = X.first.begin(); itX != X.first.end(); ++itX){
        fprintf(F,"%d ", *itX);
        //debug
        Xsize.first++;
    }

    for( std::set<node>::iterator itX = X.second.begin(); itX != X.second.end(); ++itX){
        fprintf(F,"%d ", *itX);
        //debug
        Xsize.second++;
    }
    //debug
    assert( Xsize.first == X.first.size() );
    assert( Xsize.second == X.second.size() );
    fprintf(F, "\n");
}

void print_nodeint_set(std::pair<std::vector<nodeint>,std::vector<nodeint> > *set, std::string name) {

  std::cout << name << std::endl << "\tfirst: size" << set->first.size();
    for(std::vector<nodeint>::iterator i = set->first.begin(); i != set->first.end(); ++i) {
    // for(auto i = set->first.begin(); i != set->first.end(); ++i) {
      std::cout << "(" << i->v << ",[";
      fprintf(stdout, "%g, %g]), ", i->b, i->e);
    }
    std::cout << std::endl << "\tsecond: size" << set->second.size() << " ";
    for(std::vector<nodeint>::iterator i = set->second.begin(); i != set->second.end(); ++i) {
    // for(auto i = set->first.begin(); i != set->first.end(); ++i) {
        std::cout << "(" << i->v << ",[" << i->b << "," << i->e << "]), ";
    }
    std::cout << std::endl;
}

std::vector<nodeint> inter_neighbours(std::vector<nodeint> my_nodeints, node u, float b, float e, linktype *L){
    // How to check cleverly is u is in top or bot ? Do we care ?
    std::vector<nodeint> res;

    // res = new std::vector<nodeint>();

    //loop over all (node, interval) pairs
    // Dupliquer pour second ?
    for(std::vector<nodeint>::size_type i=0; i< my_nodeints.size(); i++) {
        nodeint cur_ni = my_nodeints[i];

        //debug
        //bool show_debug = ( ((u == 891) && (cur_ni.v == 690) ) || ((u == 690) && (cur_ni.v == 635)) || ((u==635) && (cur_ni.v == 627)) || ((u==627) && (cur_ni.v == 609)) || ((u==627) && (cur_ni.v==690)) );
        //bool show_debug_fin =  ((u==627) && (cur_ni.v==690));
        //if( show_debug ){
        //fprintf(stdout, "inter_neighbours, u:  %d  [%d, %d] vs %d [%d %d]\n", u, b, e, cur_ni.v, cur_ni.b, cur_ni.e);
        //}

        std::map<node, std::vector<interval> >::iterator mymap = L[u].find(cur_ni.v);

        if( mymap != L[u].end() ){
            //cur_ni.v is a neighbour of u, need to find links that intersect with cur_ni times
            // and with [b,e]
            std::vector<interval> links = mymap->second;

            //if(show_debug_fin ){
            // for(int idebug=0; idebug < links.size(); idebug++)
            //    fprintf(stdout, "[%d %d], ", links[idebug].b, links[idebug].e);
            // fprintf(stdout, "\n");
            

            //compute intersection of [cur_ni.b, cur_ni.e] and [b,e]
            if( !( (b > cur_ni.e) || (cur_ni.b > e) ) ){
                float inter_b, inter_e;
                if( cur_ni.b > b )
                    inter_b = cur_ni.b;
                else
                    inter_b = b;
                if( cur_ni.e < e )
                    inter_e = cur_ni.e;
                else
                    inter_e = e;

                assert( inter_b <= inter_e );
		assert( (inter_e <= e) && (inter_b >= b) );
                int s = links.size();
                assert( s >= 1 );
                int first = 0;
                int last = s-1;
                int middle = (first+last) / 2;
                int found;

                while( first <= last){
                    if( links[middle].b > inter_b )
                        last = middle-1;
                    else{
                        if( links[middle].e < inter_b )
                            first = middle + 1;
                        else{
                            //found
                            assert( links[middle].b <= inter_b );
                            assert( links[middle].e >= inter_b );
			    // first link to treat after having treated links[middle], see below
                            found = middle+1;
                            break;
                        }
                    }
                    middle = (first + last)/2;
                }

                if( first > last ){
                    // inter_b is not in any link
                    assert( first == last + 1 );
                    // if last == size-1 inter_b is larger than the end of the last link
                    if( last == s-1 ){
                        assert( inter_b > links[s-1].e );
                        found = s;
                    }
                    else{
                        if( last == -1 ){
                            assert( inter_b < links[0].b );
                            found = 0;
                        }
                        else{
                            //inter_b is between two elements of links
                            assert( inter_b > links[last].b );
                            assert( inter_b < links[first].b );
                            found = first;
                        }
                    }
                }
                else{
                    //found first link that contains inter_b
                    nodeint new_ni;
                    new_ni.v = cur_ni.v;
                    new_ni.b = inter_b;
                    if( inter_e > links[middle].e )
                        new_ni.e = links[middle].e;
                    else
                        new_ni.e = inter_e;
		    assert( (new_ni.b >= b) && (new_ni.e <=e) );
                    res.push_back(new_ni);
                    //debug
                    //if( show_debug )
                    // fprintf(stdout, "  Pushing back %d [%f, %f]\n", new_ni.v, new_ni.b, new_ni.e);
                    assert( new_ni.b <= new_ni.e );
                }
                //put in res all links after the current one that intersect with [inter_b, inter_e]
                while( (found < s) && (links[found].b <= inter_e) ){
                    nodeint new_ni;
                    new_ni.v = cur_ni.v;
                    new_ni.b = links[found].b;
                    if( inter_e > links[found].e )
                        new_ni.e = links[found].e;
                    else
                        new_ni.e = inter_e;
		    assert( (new_ni.b >= b) && (new_ni.e <=e) );
                    res.push_back(new_ni);

                    //debug
                    //if( show_debug )
                    // fprintf(stdout, "  Pushing back %d [%f, %f]\n", new_ni.v, new_ni.b, new_ni.e);

                    assert( new_ni.b <= new_ni.e );

                    found ++;
                }
            }
        }
    }


    // std::vector<nodeint> res2;
    // nodeint i;
    // i.v = 1; i.b = 1; i.e = 3;
    // res2.push_back(i);
    // std::cout << "RES" << std::endl;
    // for(std::vector<nodeint>::iterator i = res.begin(); i!= res.end(); ++i)
    //     std::cout << i->v << ", [" << i->b << "," << i->e << "]" << std::endl;
    // std::cout << "RES" << std::endl;

    return res;
}

std::vector<nodeint> *inter_time(std::vector<nodeint> *my_nodeints, float b, float e){
  std::vector<nodeint> *res = new std::vector<nodeint>();

  // fprintf(stdout, "---- inter_time ------%g %g\n", b, e);
  // for(int i=0; i< my_nodeints->size();i++){
  //   fprintf(stdout, "%d [%g %g] ", my_nodeints->at(i).v, my_nodeints->at(i).b, my_nodeints->at(i).e);
  // }
  // fprintf(stdout, "\n");
  
  for(std::vector<nodeint>::size_type i=0; i<my_nodeints->size(); i++){
    nodeint cur_ni = (*my_nodeints)[i];
    if( ! (( cur_ni.e < b ) || ( e < cur_ni.b )) ){
      nodeint toinsert;
      if( cur_ni.b > b )
	toinsert.b = cur_ni.b;
      else
	toinsert.b = b;
      if( cur_ni.e < e )
	toinsert.e = cur_ni.e;
      else
	toinsert.e = e;
      toinsert.v = cur_ni.v;
      res->push_back(toinsert);
    }
  }
  // for(int i=0; i< res->size();i++){
  //   fprintf(stdout, "%d [%g %g] ", res->at(i).v, res->at(i).b, res->at(i).e);
  // }
  // fprintf(stdout, "\n");
  
  // fprintf(stdout, "------\n");
  return res;
}

//return true if and only if i is included in piv and there is a link
//between the pivot node and w during (at least) i
// This version is the same as in the non-bipartite case
bool is_in_delta_neighboorhood(nodeint piv, nodeint w_i, linktype *L){
  //  fprintf(stderr, "Entering is_in_delta_neighboorhood\n");
  if( ( piv.b <= w_i.b ) && (piv.e >= w_i.e) ){
    // i is included in pivot's interval
    // check if there is a link existing between w and pivot containing i
    std::map<node, std::vector<interval> >::iterator my_map = L[piv.v].find(w_i.v);
    if( my_map != L[piv.v].end() ){
      std::vector<interval> links = my_map->second;
      int s = links.size();
      int first = 0;
      int last = s-1;
      int middle = (first+last) / 2;
	
      while( first <= last){
	if( links[middle].b > w_i.b )
	  last = middle-1;
	else{
	  if( links[middle].e < w_i.b )
	    first = middle + 1;
	  else{
	    //found
	    assert( links[middle].b <= w_i.b );
	    assert( links[middle].e >= w_i.b );
	    break;
	  }
	}
	middle = (first + last)/2;
      }
      if( first > last ){
	//fprintf(stderr, "Leaving is_in_delta_neighboorhood\n");
	return false;
      }
      else{
	//fprintf(stderr, "Leaving is_in_delta_neighboorhood\n");
	return ( links[middle].e >= w_i.e );
      }
      
    }
    else{
      //fprintf(stderr, "Leaving is_in_delta_neighboorhood\n");
      return false;
    }
    
  }
  else{
    //fprintf(stderr, "Leaving is_in_delta_neighboorhood\n");
    return false;
  }
}

nodeint pivot(std::pair<std::vector<nodeint>, std::vector<nodeint> > *P, std::pair<std::vector<nodeint>, std::vector<nodeint> > *X, linktype *L, bool *is_in_first, bool* is_in_X){
  //fprintf(stderr, "Entering pivot\n");
  
  int max_removed = -1;
  nodeint final_pivot;
  bool piv_is_in_first, piv_is_in_X;

  //fprintf(stdout, "pivot sizes: %d %d %d %d\n", P->first.size(), P->second.size(), X->first.size(), X->second.size());
  assert( (P->first.size() > 0) || (P->second.size() > 0) );
  //loop over all nodeints in P and X
  //Pfirst
  for(std::vector<nodeint>::size_type i=0; i<P->first.size(); i++){
    nodeint piv = P->first.at(i);
    //compute number of elements that this pivot removes from P
    int nb_rem = 0;
    for(std::vector<nodeint>::size_type j=0; j<P->first.size(); j++){
      // A node should not be considered for a recursive call if it and piv are in top
      // iff this node's interval is included in piv's
      nodeint my_ni = P->first[i];
      if( (my_ni.b >= piv.b) && (my_ni.e <= piv.e) )
	// TODO the count is off by 1 in this case because piv is considered as removed here
	nb_rem++;
    }
    for(std::vector<nodeint>::size_type j=0; j<P->second.size(); j++){
      if( is_in_delta_neighboorhood(piv, P->second.at(j), L) )
	nb_rem++;
    }
    // Number of nodes removed are those removed in P->second, and all nodes in P->first except pivot
    assert( nb_rem + (int)P->first.size() -1 >= 0 );
    if( nb_rem > max_removed ){
      max_removed = nb_rem;
      //fprintf(stdout, "max removed %d\n%",max_removed);
      final_pivot = piv;
      piv_is_in_first = true;
      piv_is_in_X = false;
    }
    
  }

  //Psecond
  for(std::vector<nodeint>::size_type i=0; i<P->second.size(); i++){
    nodeint piv = P->second[i];
    int nb_rem = 0;
    for(std::vector<nodeint>::size_type j=0; j<P->first.size(); j++){
      if( is_in_delta_neighboorhood(piv, P->first[j], L) )
	nb_rem++;
    }
    for(std::vector<nodeint>::size_type j=0; j<P->second.size(); j++){
      // A node should not be considered for a recursive call if it and piv are in second
      // iff this node's interval is included in piv's
      nodeint my_ni = P->second[i];
      if( (my_ni.b >= piv.b) && (my_ni.e <= piv.e) )
	// TODO the count is off by 1 in this case because piv is considered as removed here
	nb_rem++;
    }
    assert( nb_rem >= 0 );
    if( nb_rem  > max_removed ){
      max_removed = nb_rem;
      final_pivot = piv;
      piv_is_in_first = false;
      piv_is_in_X = false;
    }
  }
  
  //Xfirst
  for(std::vector<nodeint>::size_type i=0; i<X->first.size(); i++){
    nodeint piv = X->first.at(i);
    //compute number of elements that this pivot removes from P
    int nb_rem = 0;
    for(std::vector<nodeint>::size_type j=0; j<P->first.size(); j++){
      // A node should not be considered for a recursive call if it and piv are in top
      // iff this node's interval is included in piv's
      nodeint my_ni = P->first[i];
      if( (my_ni.b >= piv.b) && (my_ni.e <= piv.e) )
	// TODO the count is off by 1 in this case because piv is considered as removed here
	nb_rem++;
    }
    for(std::vector<nodeint>::size_type j=0; j<P->second.size(); j++){
      if( is_in_delta_neighboorhood(piv, P->second.at(j), L) )
	nb_rem++;
    }
    if( nb_rem > max_removed ){
      max_removed = nb_rem;
      final_pivot = piv;
      piv_is_in_first = true;
      piv_is_in_X = true;
    }
  }


  //Xsecond
  for(std::vector<nodeint>::size_type i=0; i<X->second.size(); i++){
    nodeint piv = X->second.at(i);
    //compute number of elements that this pivot removes from P
    int nb_rem = 0;
    for(std::vector<nodeint>::size_type j=0; j<P->first.size(); j++){
      if( is_in_delta_neighboorhood(piv, P->first.at(j), L) ){
	nb_rem++;
      }
    }
    for(std::vector<nodeint>::size_type j=0; j<P->second.size(); j++){
      // A node should not be considered for a recursive call if it and piv are in second
      // iff this node's interval is included in piv's
      nodeint my_ni = P->second[i];
      if( (my_ni.b >= piv.b) && (my_ni.e <= piv.e) )
	// TODO the count is off by 1 in this case because piv is considered as removed here
	nb_rem++;
    }

    if( nb_rem  > max_removed ){
      max_removed = nb_rem;
      final_pivot = piv;
      piv_is_in_first = false;
      piv_is_in_X = true;
    }
  }
  
  assert( max_removed != -1 );
  //fprintf(stderr, "Leaving pivot\n");
  //fprintf(stdout, "Nb removed %d\n", max_removed);
  *is_in_first = piv_is_in_first;
  *is_in_X = piv_is_in_X;
  return final_pivot;
}


  //depth is used for debug
void bron_kerbosch(std::pair<std::vector<nodeint>, std::vector<nodeint> > *P, std::pair<std::vector<nodeint>, std::vector<nodeint> > *Pmax, Clique R, std::pair<std::vector<nodeint>, std::vector<nodeint> > *X, std::pair< std::vector<nodeint>, std::vector<nodeint> > *Xmax, linktype *L, int depth){
        
  std::pair< std::vector<nodeint>, std::vector<nodeint> > *Pprimemax;
  std::pair< std::vector<nodeint>, std::vector<nodeint> > *Xprimemax;
  //TODO: use other structure than vector for P, to allow efficient deletion?

  //if( (Pmax->size() == 0) && (Xmax->size() == 0) ){
  if( ((Pmax->first.size() == 0) && (Xmax->first.size() == 0 )) && ((Pmax->second.size() == 0) && (Xmax->second.size() == 0)) ) {
    if( R.X.first.size() >= 1 || R.X.second.size() >= 1) {
      // if( true ){
      // std::cout << "is maximal :" <<std::endl;
      R.print_clique(stdout);
      fflush(stdout);
      // Res->insert(R);
    }
  }

  unsigned int ifirst = 0;
  unsigned int isecond = 0;
  nodeint piv;
  bool piv_is_in_first, piv_is_in_X;
  if( ( ifirst < P->first.size() ) || ( isecond < P->second.size() ) ){
    assert( (P->first.size() > 0) || (P->second.size() > 0) );
    //fprintf(stdout, "Pivoting in %d %d %d %d\n", ifirst, P->first.size(), isecond, P->second.size() );
    if( depth <=3){
      fprintf(stdout, "---- BK, depth = %d \n", depth);
      fprintf(stdout, "Current clique: ");
      R.print_clique(stdout);
      print_nodeint_set(P, "P");
      print_nodeint_set(Pmax, "P_max");
      print_nodeint_set(X, "X");
      print_nodeint_set(Xmax, "X_max");

    }


    piv = pivot(P, X, L, &piv_is_in_first, &piv_is_in_X);

    if( depth <= 3){
      fprintf(stdout, "Pivot : %d [%g %g]\n", piv.v, piv.b, piv.e);
    }
    // fprintf(stdout, "Current clique: ");
    // R.print_clique(stdout);
    // fprintf(stdout, "Pivot: %d [%g %g]\n", piv.v, piv.b, piv.e);
    // if( (piv.v > 20) || (piv.v <0) ){
    //   fprintf(stdout, "wrong pivot\n");
    //   print_nodeint_set(P, "P ");
    // }

  }

  //  while( P->first.size() != 0 || P->second.size() != 0){
  while( (ifirst < P->first.size()) || (isecond < P->second.size()) ){
    //  for(int iP=0; iP<P->size(); iP ++){
    nodeint cur_ni;
    bool side; // 0 for top, 1 for bot
    
    if( ifirst < P->first.size()) {
      cur_ni = P->first.at(ifirst);
      side = 0;
    }
    else {
      cur_ni = P->second.at(isecond);
      side = 1;
    }
    assert(side == 0 || side == 1);

    bool recursive_call = true;
    // TODO this should be double-checked
    // nodeints for which a recursive calls should be made if piv is in first are:
    // - piv if piv is in P (and not in X)
    // - elements of P->first for which the time interval is not included in piv's interval -> IMPLEMENT THIS BELOW
    // - elements of P->second that are not in the neighboorhood of P
    // The order between the two ifs below is important: piv should be considered
    // even though it is in the same side as piv

    if( piv_is_in_first ){
      if( side == 0 ){
	if( (cur_ni.b >= piv.b) && (cur_ni.e <= piv.e) )
	  recursive_call = false;
	else
	  recursive_call = true;
      }
      else{
	recursive_call = !is_in_delta_neighboorhood(piv, cur_ni, L);
      }
    }
    else{
      if( side == 1 )
	//piv is in second. Perform a recursive call only if cur_ni's interval is not included in piv's
	if( (cur_ni.b >= piv.b) && (cur_ni.e <= piv.e) )
	  recursive_call = false;
	else
	  recursive_call = true;
      else
	recursive_call = !is_in_delta_neighboorhood(piv, cur_ni, L);
    }
    if( ! piv_is_in_X ){
      if( (cur_ni.v == piv.v) && (cur_ni.b == piv.b) && (cur_ni.e == piv.e) )
	//if piv is in P and piv == cur_ni -> recursive_call (regardless of side)
	recursive_call = true;
    }

    if( recursive_call ){
      if(side == 0) {
	P->first.erase(P->first.begin()+ifirst);
	//X->first.push_back(cur_ni);
      }
      else {
	P->second.erase(P->second.begin()+isecond);
	//X->second.push_back(cur_ni);   
      }
      
      
      if( depth <=3){
	fprintf(stdout, "Taking (node, interval) pair: %d [%f, %f]\n", cur_ni.v, cur_ni.b, cur_ni.e);
      }
      Clique Rprime = Clique(R);
      
      if(side == 0) {
	Rprime.X.first.insert( cur_ni.v );
      }
      else {
	Rprime.X.second.insert( cur_ni.v );
      }
      
      Rprime.b = cur_ni.b;
      Rprime.e = cur_ni.e;
      std::pair< std::vector<nodeint>, std::vector<nodeint> > *Pprime = new std::pair< std::vector<nodeint>, std::vector<nodeint> >();
      std::pair< std::vector<nodeint>, std::vector<nodeint> > *Xprime = new std::pair< std::vector<nodeint>, std::vector<nodeint> >();
      
      if(side == 0) {                
	Pprime->second = inter_neighbours(P->second, cur_ni.v, cur_ni.b, cur_ni.e, L); // for cur_ni in top, returns the (u,I) in bot in cur_ni's neighborhood.
	// TODO Ugly
	std::vector<nodeint> *res = inter_time(&(P->first), cur_ni.b, cur_ni.e);
	Pprime->first = *res;
	delete res;
	Xprime->second = inter_neighbours(X->second, cur_ni.v, cur_ni.b, cur_ni.e, L);
	res = inter_time(&(X->first), cur_ni.b, cur_ni.e);
	Xprime->first = *res;
	delete res;
      }
      else {
	Pprime->first = inter_neighbours(P->first, cur_ni.v, cur_ni.b, cur_ni.e, L);
	std::vector<nodeint> *res = inter_time(&(P->second), cur_ni.b, cur_ni.e);
	Pprime->second = *res;
	delete res;
	Xprime->first = inter_neighbours(X->first, cur_ni.v, cur_ni.b, cur_ni.e, L);
	res = inter_time(&(X->second), cur_ni.b, cur_ni.e);
	Xprime->second = *res;
	delete res;
	
      }
      
      //assert( (Pprime->first.size() < P->first.size()) || (Pprime->second.size() < P->second.size() ) );
      
      
      //TODO: Compute Pprimemax and Xprimemax on the fly, at the same time as Pprime and Xprime
      Pprimemax = new std::pair< std::vector<nodeint>, std::vector<nodeint> >;
      
      // TODO: Factorize duplicate code
      for(std::vector<nodeint>::size_type iprime=0; iprime < Pprime->second.size(); iprime ++){
	// TOD : commented lines not relevant?
	//std::cout << "debug: cur_ni, v=" << cur_ni.v << ", b=" << cur_ni.b << ", e=" << cur_ni.e << std::endl; 
	//std::cout << "debug: P', v=" << Pprime->first[iprime].v << ", b=" << Pprime->first[iprime].b << ", e=" << Pprime->first[iprime].e << std::endl; 
	if( ( (Pprime->second)[iprime].b == cur_ni.b ) && ( (Pprime->second)[iprime].e == cur_ni.e ) )
	  Pprimemax->second.push_back( (Pprime->second)[iprime] );
      }
    
    
      for(std::vector<nodeint>::size_type iprime=0; iprime < Pprime->first.size(); iprime ++){
	if( ( (Pprime->first)[iprime].b == cur_ni.b ) && ( (Pprime->first)[iprime].e == cur_ni.e ) )
	  Pprimemax->first.push_back( (Pprime->first)[iprime] );
      }
      
      
      Xprimemax = new std::pair< std::vector<nodeint>, std::vector<nodeint> >;
      for(std::vector<nodeint>::size_type iprime=0; iprime < Xprime->first.size(); iprime ++){
	if( ( (Xprime->first)[iprime].b == cur_ni.b ) && ( (Xprime->first)[iprime].e == cur_ni.e ) )
	  Xprimemax->first.push_back( (Xprime->first)[iprime] );
      }
      for(std::vector<nodeint>::size_type iprime=0; iprime < Xprime->second.size(); iprime ++){
	if( ( (Xprime->second)[iprime].b == cur_ni.b ) && ( (Xprime->second)[iprime].e == cur_ni.e ) )
	  Xprimemax->second.push_back( (Xprime->second)[iprime] );
      }
      // if(depth <= 2 ){
      // 	std::cout << "Calling BK, depth " << depth <<  " with :" <<std::endl << "    ";
      // 	Rprime.print_clique(stdout);
      // 	print_nodeint_set(Pprime, "P'");
      // 	print_nodeint_set(Pprimemax, "P'_max");
      // 	print_nodeint_set(Xprime, "X'");
      // 	print_nodeint_set(Xprimemax, "X'_max");
      // }
      //char cont = 'd';
      // while(cont != 'c') {
      //     std::cin.get();
      //     std::cin >> cont;
      //     if(cont == 'c') break;
      //     std::cout << "cont is " << cont << std::endl;
      // }
      bron_kerbosch(Pprime, Pprimemax, Rprime, Xprime, Xprimemax, L, depth+1);
      delete Pprime;
      delete Xprime;
      delete Pprimemax;
      delete Xprimemax;
      if(side == 0) {
	X->first.push_back(cur_ni);
      }
      else {
	X->second.push_back(cur_ni);   
      }
    }
    else{
      //increment ifirst or isecond
      if( side == 0) {
	ifirst++;
      }
      else {
	isecond++;
      }

    }
  }
}

//Reads the file and populates V with all nodes seen
std::pair< std::set<node>, std::set<node> > initialize(std::pair< std::set<node>, std::set<node> > *V, linktype *L, int n, float *alpha, float *omega){
        char line[MAX_LINE_LENGTH];
        float b,e;
        node u, v;

        FILE* f=stdin;

        fgets(line,MAX_LINE_LENGTH,f);

        bool first = true;
        while ( !feof(f) ){
            sscanf(line, "%f %f %d %d\n", &b, &e, &u, &v);
            //update alpha
            if( first ){
                *alpha = b;
                first = false;
            }
            //update omega
            if( e > *omega )
                *omega = e;

            //fprintf(stdout, "%s", line);
            //update link structure
            // if node is not in V we need to initialize its map
            assert( u < n );
            assert( v < n );
            bool is_in = V->first.find(u) != V->first.end();
            if( !is_in ){
                L[u] = linktype();
            }
            is_in = V->second.find(v) != V->second.end();
            if( !is_in ){
                L[v] = linktype();
            }

            interval inter;
            inter.b = b;
            inter.e = e;

            //Add the current link in u's map
            std::map<node, std::vector<interval> >::iterator mymap = L[u].find(v);

            if( mymap == L[u].end() ){
                //need to insert v in u's map with a vector containing just the single current interval
                L[u].insert(std::pair<node,std::vector<interval> >(v, std::vector<interval>(1,inter)));
            }
            else{
                mymap->second.push_back(inter);
            }

            //Add the current link in v's map
            mymap = L[v].find(u);
            if( mymap == L[v].end() ){
                L[v].insert(std::pair<node, std::vector<interval> >(u, std::vector<interval>(1,inter)));
            }
            else{
                mymap->second.push_back(inter);
            }

            //update V
            V->first.insert(u);
            V->second.insert(v);
            fgets(line,MAX_LINE_LENGTH,f);
        }

        return *V;
    }

    //For debug
    void print_linkstream(std::vector<linktype> L){
        for(std::vector<nodeint>::size_type i=0; i< L.size(); i++){
            fprintf(stdout, "%d : ", i);
            for(std::map<node, std::vector<interval> >::iterator it = L[i].begin(); it!=L[i].end(); ++it){
                node neigh = it->first;
                for(std::vector<nodeint>::size_type iu=0; iu < it->second.size(); iu++){
                    fprintf(stdout, "%d (%f %f) ", neigh, it->second[iu].b, it->second[iu].e);
                }
                fprintf(stdout, "\n");
            }
        }
    }

    //For debug
    void print_cliqueset(std::set<Clique> mySet){
        int setsize = 0;

        std::set<Clique>::iterator it = mySet.begin();
        Clique cur_c = *it;
        Clique prev_c = *it;
        while( it != mySet.end() ){
            //  for( std::set<Clique>::iterator it = mySet.begin(); it != mySet.end(); ++it ){
            it->print_clique(stdout);
            setsize ++;
            ++it;
            cur_c = *it;
            //fprintf(stderr, "print_cliqueset:\n prev_c: ");
            //prev_c.print_clique(stderr);
            //fprintf(stderr, "cur_c: ");
            //cur_c.print_clique(stderr);
            assert ( (it == mySet.end()) || (prev_c < cur_c) );
            //assert( prev_c.b <= cur_c.b );
            prev_c = *it;
        }

        assert( setsize == mySet.size() );
        }

int main(int argc, char **argv){
  // By convention, in pairs of vectors of nodes, first is top, second is bot. 
  std::pair< std::set<node>, std::set<node> > V;
  //  std::vector<nodeint> *Pinit;
  std::pair< std::vector<nodeint>, std::vector<nodeint> > Pinit;
  // nodeint i;
  // i.v = 0;
  // i.b = 2;
  // i.e = 3;
  // Pinit.first.push_back(i);
  
  // for( std::vector<nodeint>::iterator i = Pinit.first.begin(); i != Pinit.first.end(); ++i)
  //     std::cout << i->v << " [" << i->b << "," << i->e << "]" << std::endl;
  
  // Pinit.first = std::vector<nodeint>();
  // Pinit.second = std::vector<nodeint>();
  
  std::pair< std::vector<nodeint>, std::vector<nodeint> > Pmax;
  std::pair< std::vector<nodeint>, std::vector<nodeint> > Xmax;
  std::pair< std::vector<nodeint>, std::vector<nodeint> > X;
  linktype* newLinks;
  int n;
  float alpha, omega;
  bool debug = false;
  int depth = 100;
  
  if( argc != 6 ){
    fprintf(stderr, "Usage: %s n b e u v\n", argv[0]);
    return(1);
  }
  
  sscanf(argv[1], "%d", &n);
  //assert( (newLinks = (linktype*)malloc(n*sizeof(linktype))) != NULL );
  //for(int i=0;i<n;i++)
  //newLinks[i] = linktype();
  newLinks = new linktype[n];
  
  alpha = 0.0;
  omega = 0.0;
  // TODO Keep track of actual value of n in initialize
  V = initialize(&V, newLinks, n, &alpha, &omega);
  
  //fprintf(stdout, "Initialize done\n");
  //Compute initial P from newLinks
  //  *Pinit = std::vector<nodeint>();
  
  for(int v=0; v<n; v++){
    
    nodeint cur_ni;
    cur_ni.v = v;
    cur_ni.b = alpha;
    cur_ni.e = omega;
    
    bool is_in_top = V.first.find(v) != V.first.end();
    if( is_in_top ) {
      Pinit.first.push_back(cur_ni);
    }
    else {
      Pinit.second.push_back(cur_ni);
    }
  }
  
  //find first link in dataset
  //find first node in first
  // std::set<node>::iterator V_it = V.first.begin();
  // node u = *V_it;
  node u;
  sscanf(argv[4], "%d", &u);
  node v;
  sscanf(argv[5], "%d", &v);
  //bool is_in_top = V.first.find(u) != V.first.end();
  //assert( is_in_top );
  // find first neighbour of u
  //std::map<node, std::vector<interval> >::iterator it = newLinks[u].begin();
  //assert( it != newLinks[u].end() );
  //nodeint u_first_ni;
  //u_first_ni.v = it->first;
  //u_first_ni.b = it->second[0].b;
  //u_first_ni.e = it->second[0].e;
  float b;
  sscanf(argv[2], "%g", &b);
  float e;
  sscanf(argv[3], "%g", &e);
  
  //fprintf(stderr, "P done\n");
  
  //TODO Assumption here that no links lasts as long as [alpha, omega]
  //  *Pmax=std::vector<nodeint>();
  //  *X=std::vector<nodeint>();
  //*Xmax=std::vector<nodeint>();
  Clique R = Clique(alpha, omega, std::pair<std::set<node>, std::set<node> >());
  
  // Add first link to R
  R.X.first.insert(u);
  R.X.second.insert(v);
  R.b = b;
  R.e = e;
  
  //R.print_clique(stdout);
  
  std::vector<nodeint> Pdeuxs = inter_neighbours(Pinit.second, u, R.b, R.e, newLinks);
  std::vector<nodeint> Pdeuxf = inter_neighbours(Pinit.first, v, R.b, R.e, newLinks);
  //Need to remove u and v from Pdeux
  for(std::vector<nodeint>::size_type i=0; i< Pdeuxf.size(); i++){
    nodeint cur_ni = Pdeuxf[i];
    if( (cur_ni.v == u) && (cur_ni.b == R.b) && (cur_ni.e == R.e) )
      Pdeuxf.erase(Pdeuxf.begin()+i);
  }
  for(std::vector<nodeint>::size_type i=0; i< Pdeuxs.size(); i++){
    nodeint cur_ni = Pdeuxs[i];
    if( (cur_ni.v == v) && (cur_ni.b == R.b) && (cur_ni.e == R.e) )
      Pdeuxs.erase(Pdeuxs.begin()+i);
  }
  
  std::pair<std::vector<nodeint>, std::vector<nodeint> >Pdeux = std::pair<std::vector<nodeint>, std::vector<nodeint> >(Pdeuxf, Pdeuxs);
  //Pmax for start
  for(std::vector<nodeint>::size_type i=0; i < Pdeux.first.size(); i++){
    if( ( (Pdeux.first)[i].b == R.b ) && ( (Pdeux.first)[i].e == R.e ) )
      Pmax.first.push_back( (Pdeux.first)[i] );
  }
  for(std::vector<nodeint>::size_type i=0; i < Pdeux.second.size(); i++){
    if( ( (Pdeux.second)[i].b == R.b ) && ( (Pdeux.second)[i].e == R.e ) )
      Pmax.second.push_back( (Pdeux.second)[i] );
  }
  
  
  // debug
  if( debug ){
    depth = -100;
    fprintf(stdout, "BK début, lien %g %g %d %d\n", R.b, R.e, u, v);
    print_nodeint_set(&Pdeux, "Pdeux");
    print_nodeint_set(&Pmax, "Pmax");
    print_nodeint_set(&X, "X");
    print_nodeint_set(&Xmax, "Xmax");
  }
  
  bron_kerbosch(&Pdeux, &Pmax, R, &X, &Xmax, newLinks, depth);
}
