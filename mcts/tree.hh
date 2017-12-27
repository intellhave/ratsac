#ifndef TREE_HH
#define TREE_HH

#include <vector>
#include <math.h>
#include "State.hh"
#include "node.hh"


namespace MCTS{
  
  class Tree{
  private:
    int totalNodes;
    

  public:
    Node *root;
    Tree(state* _state);
    ~Tree();
    void addChild(Node *parent, Node *child);
    Node* expand(Node *node);
    Node* bestChild(Node *node, double Cp);
    void backPropagate(Node *node, int reward);            
  };
}



#endif
