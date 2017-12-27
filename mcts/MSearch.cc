#include "MSearch.hh"

using namespace std;
using namespace Eigen;

namespace MCTS{
    
  MSearch::MSearch(MatrixXd *x_, VectorXd *y_, const int &nP_, const double &th_){
    xData = x_;
    yData = y_;
    th = th_;

    N = xData->rows();
    nPoints = nP_;
    if (nPoints == N)
      problemType = LINEAR;
    else if (nPoints*4 == N)
      problemType = QUASICONVEX;
    else {      
      assert(false);
    }
    dim = xData->cols();
    treeTravelIter = 1;
    
  }


  void MSearch::startSearch(){
    
    VectorXd theta0f = Eigen::VectorXd::Random(dim); //Dummy initialization
    
    set<int>vSet; 
    set<int>sSet;

    vSet.clear(); //Violation Set
    
    for (int i=0; i<nPoints; i++) //Get all the points to support set;
      sSet.insert(i);
    //Initialize static variables
    state::stateMap.clear();
    state::bestTheta.resize(dim);
    state::problemType = problemType;
    state::AA = AA;
    state::bb = bb;
    state::cc = cc;
    state::dd = dd;
 
    clock_t tick = clock();
    state initState(xData, yData, vSet, sSet, th, MCTS::DOWN, theta0f);
    
    Tree *searchTree = new Tree(&initState); //Create a tree with the initial state
    state::totalIter = 0;
    state::maxInls = 0;
    
    for (int treeit = 0; treeit < treeTravelIter; treeit++){      
      Node *currentNode = searchTree->root; //Start from the root node
      
      while (!currentNode->nodeState->terminated()){ //While the terminated state is not reached
	state *currentState = currentNode->nodeState;
	int maxChildren = currentState->basics.size();
	Node *prevNode = currentNode;
	int iter = 9;
	currentNode = uctSearch(searchTree, currentNode, iter);
	
	if (prevNode == currentNode){
	  cout << "---------=======-------------"<<endl;
	  cout << currentNode->nodeState->ptsInls << "\n";	  
	  break;
	}	
      }

      theta = state::bestTheta;
      cout << currentNode->nodeState->minimax <<"\n";
      cout << currentNode->nodeState->supportSet.size()<<"\n";
           
    }
  }
  
}

