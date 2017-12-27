#include <iostream>
#include "mcts.hh"

namespace MCTS{
  Node* uctSearch(Tree* searchTree, Node* node, int iter){
    node->visits = 0;
    for (int i=0; i<iter; i++){
      Node *expandedNode = treePolicy(searchTree, node);
      int reward = defaultPolicy(expandedNode->nodeState);
      searchTree->backPropagate(expandedNode, reward);
    }
    return searchTree->bestChild(node,0);
  }

  Node* treePolicy(Tree *searchTree, Node* node){
    while (!node->nodeState->feasible()){
      if (!node->fullyExpanded()){
	return searchTree->expand(node);
      }
      else{
	//double Cp = 1/sqrt(2);
	double Cp = 500;
	Node *prevNode = node;
	node=searchTree->bestChild(node, Cp);
	if (prevNode == node)
	  return node;
      }	
    }
    return node;
  }
  
  int defaultPolicy(state *state){    
    MCTS::state *cState = state->nextRandomState();
    int rewd = 0;    
    while (!cState->terminated()){
      MCTS::state* pState = cState;
      cState = cState->nextRandomState();
      int sInls = cState->inls;
      if (sInls  > rewd && sInls < 100000)
	rewd = sInls;      
      /*      if (pState!=state && pState!=cState)
	      delete pState;
	      else
	      break;
      */
    }    
    return rewd;
  }
  
}
  
