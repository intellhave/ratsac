

#include "node.hh"

namespace MCTS{

  Node::Node(state *_state){
    visits = 1;
    reward = 0;
    nChildren = 0;
    expandedAll = false;
    children.clear();
    nodeState = _state;    
  }
  
  Node::~Node(){
    //    delete nodeState;
  }

  void Node::updateReward(int _reward){
    //reward =max(reward,_reward);
    reward += _reward;
  }

  bool Node::fullyExpanded(){
    if (  (nodeState->direction == DOWN && children.size()==nodeState->basics.size())
	  ||(nodeState->direction == UP && children.size()==nodeState->violationSet.size())  )
      return true;
    else
      return false;
    
  }
  
}
