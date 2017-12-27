#ifndef NODE_HH
#define NODE_HH

#include <vector>
#include "State.hh"


namespace MCTS{
  
  class Node{
  public:
    int visits;
    int reward;
    int  nChildren;
    bool expandedAll;
    int bestInlierSize;
    std::vector<Node*> children;
    Node* parent;
    state *nodeState; //Each node corresponds to a state;

    Node(state *_state);
    ~Node();
    void updateReward(int _reward);
    bool fullyExpanded();
  };

}
#endif
