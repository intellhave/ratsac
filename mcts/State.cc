#include <cstdlib>
#include <iostream>
#include "State.hh"
#include "LinearFunctions.hh"

using namespace std;
//using namespace arma;

namespace MCTS{

  //Declaration of the static variables
  int state::totalIter;
  unordered_map<std::string, state*> state::stateMap;
  int state::maxInls;
  Eigen::VectorXd state::bestTheta;
  pType state::problemType;

  Eigen::MatrixXd *state::AA;
  Eigen::MatrixXd *state::bb;
  Eigen::MatrixXd *state::cc;
  Eigen::VectorXd *state::dd;
  
  //Constructor 
  state::state(Eigen::MatrixXd *xData_, Eigen::VectorXd *yData_, const set<int>&vSet, const set<int> &sSet, double th_, stateDirection prevDirection, const Eigen::VectorXd &theta0){

    totalIter = totalIter + 1;
    
    xData = xData_;
    yData = yData_;

    N = xData->rows(); /*Data size*/   
    d = xData->cols(); /*Data dimension*/
    th = th_;
    int nc = 1;        /*Number of linear constraints per data instance;*/
    bool abs = false;
    if (problemType == LINEAR){
      state_th = th;
      abs = true;
      nc = 1;
      nPoints = N;
    }
    else{
      state_th = 0.001;
      abs = false;
      nc = 4;
      nPoints = N/4;
    }
    violationSet = vSet;
    supportSet = sSet;
    
    /*From the coverage set, extract the set of linear constraints*/
    constraintSet.clear();
    if (problemType == LINEAR)  {
      constraintSet = sSet;   /*Constraint set is the support set for linear case*/
    }
    else {
      /*Else, for each point, there is 4 equivalent constraints*/
      for (auto it = sSet.begin(); it!=sSet.end(); it++){
     	for (int rIdx = 4*(*it); rIdx < 4*(*it)+4; rIdx++)
     	  constraintSet.insert(rIdx);	
      }
    }
    /*Extract contraints to a new matrix*/
    //constraintSet = sSet;
    Eigen::MatrixXd xFit(constraintSet.size(),d);
    Eigen::VectorXd yFit(constraintSet.size());
    set<int> dataBasics;
    vector<int> idx;
    int i=0;    
    for (set<int>::iterator it = constraintSet.begin(); it!=constraintSet.end(); it++){
      idx.push_back(*it);
      xFit.row(i) = (*xData).row(*it);
      yFit[i] = (*yData)(*it);
      i++;	  
    }
    
    /*Perform fitting*/
    if (problemType==LINEAR) abs = true; else abs = false;
    xn = Chebyshev::Descend(xFit, yFit,  th, theta0, dataBasics, minimax, abs);
        
    //Extract basic set
    basics.clear();
    for (set<int>::iterator it = dataBasics.begin(); it!=dataBasics.end(); it++){      
      basics.insert(idx[*it]/nc) ;
    }
   
    /*Compute number of inliers for the current fit*/
    if (problemType==LINEAR){
      inls = findInliers(*xData, *yData, xn, state_th, true);
    }
    else{
      inls = QuasiconvexTools::findInliers(*AA, *bb, *cc, *dd, xn, th);
    }
    
    if (inls > maxInls){
      maxInls = inls;
      bestTheta = xn;
      cout << "Max Inls = " << maxInls << endl;
      //test
      // int pointsInls; 
      // if (problemType == QUASICONVEX){
      // 	ptsInls = QuasiconvexTools::findInliers(*AA, *bb, *cc, *dd, xn, th);
      // 	cout << "P Inls = " << ptsInls <<"\n";
      // }      
    }    
    direction = DOWN;
    if (prevDirection==UP){
      if (minimax <= state_th)
	direction = UP;
    }		    
  }
    
  //Check if a state is terminated (can't be descent down)
  bool state::terminated(){
    int sSize = supportSet.size();
    if (
	(minimax <= state_th) ||
	(sSize<=d+1) ||
	(sSize < maxInls)
	)
      return true;
    else
      return false;    
  }

  bool state::feasible(){
    if (((minimax<=state_th) ) || supportSet.size()<=d+1)
      return true;
    else
      return false;
  }
  
  //Generate a random state
  state* state::nextRandomState(){
    if (direction == DOWN)
      return nextRandomStateDown();	
    else
      return nextRandomStateUp();	    	
  }

  //Generate a 
  state* state::nextUnexploredState() {
    if (direction == DOWN)	    
      return nextUnexploredStateDown();
    else
      return nextUnexploredStateUp();	
  }
  
  state* state::nextRandomStateDown(){
    set<int> newV  = violationSet;
    set<int> newS = supportSet;
    
    int pickedBasic = rand() % basics.size();
    set<int>::iterator it = basics.begin();
    std::advance(it, pickedBasic);
    //Remove basic from the support set and add to violation set
    newV.insert(*it);
    set<int>::iterator erasePos = newS.find(*it);
    if (erasePos != newS.end())
      newS.erase(erasePos);

    //Check if violation is already in the map
    string vSetString = stringFromSet(newV);
    unordered_map<string, state*>::iterator mapPos = stateMap.find(vSetString);
    state *newState;
    if (mapPos != stateMap.end()){
      newState = mapPos->second;
    }
    else{      
      newState = new state(xData, yData, newV, newS, th, DOWN, xn);
      stateMap.insert(std::pair<string, state*>(vSetString, newState));
    } 
    
    return newState;
  } 
  
  state* state::nextRandomStateUp(){
	
  }
  
  state *state::nextUnexploredStateDown(){
    set<int>::iterator it;	
    for (it = basics.begin(); it!=basics.end(); it++){
      if (basicExplored.find(*it)==basicExplored.end()){		
	basicExplored.insert(*it);		
	break;
      }
    }
    //Now, remove the basic and generate a new state down
    //set<int> newV = violationSet;
    set<int> newV;
    newV = violationSet;
    set<int> newS = supportSet;
    newV.insert(*it);
    set<int>::iterator erasePos = newS.find(*it);
    if (erasePos != newS.end())
      newS.erase(erasePos);

    string vSetString = stringFromSet(newV);
    unordered_map<string, state*>::iterator mapPos = stateMap.find(vSetString);
    state *newState;
    
    if (mapPos != stateMap.end()){
      newState = mapPos->second;    
    }
    else{      
      newState = new state(xData, yData, newV, newS, th, DOWN, xn);
      stateMap.insert(std::pair<string, state*>(vSetString, newState));      
    }

    
    //Now, check if the state hash by the violation set is is the stateMap;    
    return newState;
    //return new state(xData, yData, newV, newS, th, DOWN, xn);	    	    	
  }
    
  state* state::nextUnexploredStateUp(){
    
  }

  int state::computeReward(){
    int rewd;
    if (minimax <= state_th)
      rewd =  nPoints - violationSet.size();
    else
      rewd = 0;

    if (rewd < 0 || rewd > 1000000){
      rewd = 0;
    }

    return rewd;
  }


  string stringFromSet(set<int> vSet){
    
    string str = "";
    for (set<int>::iterator it = vSet.begin(); it!=vSet.end(); it++)
      str = str+std::to_string(*it);
    return str;
    
  }
  
    
}






