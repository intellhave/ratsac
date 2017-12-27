
#include "tree.hh"

namespace MCTS{

  Tree::Tree(state *_state){
    
    root = new Node(_state);
    root->parent = NULL;
    
  }

  Tree::~Tree(){
    delete root;
  }

  void Tree::addChild(Node *parent, Node *child){
    parent->nChildren++;
    parent->children.push_back(child);
    child->parent = parent;    
  }
  
  Node* Tree::expand(Node *node){
    state* currentState = node->nodeState;
    state* nextState = currentState->nextUnexploredState();    
    Node* childNode = new Node(nextState);
    addChild(node, childNode);
    /*After adding a child, check if all the children are added*/
    if (node->nChildren == node->nodeState->basics.size())
      node->expandedAll = true; //This flags prevent more child to be added
    return childNode;
  }
  
  Node* Tree::bestChild(Node *node, double Cp){
    Node *bestChild = node;
    double bestScore = 0;
    for (int i=0; i<node->children.size(); i++){
      Node* child = node->children[i];

      if (child == NULL) continue;
      /*If the consensus size of the node is less 
	than best results, ignore it*/      
      if (child->nodeState->terminated()){
	//delete child;
	//node->children[i] = NULL;
	continue;
      }
      
      double exploit = (child->reward*1.0)/(child->visits*1.0);
      //double explore = sqrt(2*log(node->visits)/child->visits);
      double explore = sqrt(2*sqrt(node->visits)/child->visits);
      double score = exploit + Cp*explore;
      //      cout << i << "\t" << child<< "  exploit  = " << exploit << "  explore =  "  << explore
      //<< "\t" << node->visits <<"\t" << child->visits<< "\n";
      if (score > bestScore){
	bestScore = score;
	bestChild = child;	
      }
    }
    
    //If the support set of the next state is less than current best inls set
    //if (bestChild->nodeState->supportSet.size()<=state::maxInls)
    //  bestChild = node;
    return bestChild;
  }
  
  void Tree::backPropagate(Node *node, int reward){    
    node->visits++;
    node->updateReward(reward);
    while (node->parent != NULL){
      node = node->parent;
      node->visits++;
      node->updateReward(reward);  
    }
  }
  
}
