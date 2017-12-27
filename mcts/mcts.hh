#ifndef MCTS_HH
#define MCTS_HH


#include "tree.hh"
#include "State.hh"
#include "node.hh"

namespace MCTS{  
  Node* uctSearch(Tree* searchTree, Node* node, int iter);
  Node* treePolicy(Tree* searchTree, Node* node);
  int defaultPolicy(state *state);
  //void mctsSearch(const mat &x, const vec &y, const double th);

  
  
}



#endif
