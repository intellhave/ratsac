#ifndef STATE_HH
#define STATE_HH

#include <cstdlib>
#include <vector>
#include <set>
//#include <armadillo>

#include "extra.hh"
#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "quasiconvexFunctions.hh"

#include <Eigen/Dense>
#include <Eigen/Core>

#include <string>
#include <unordered_map>

using namespace std;
//using namespace arma;
using namespace LinearTools;
using namespace ToolBox;

namespace MCTS{
  enum stateDirection {UP, DOWN};
  enum pType {LINEAR, QUASICONVEX}; //Problem Type
  //Functions for hasing
  std::string stringFromSet(set<int> vSet);

  class state{

  private:
    
    set<int> basicExplored;
    int N, nPoints, d;
    //mat *x; vec *y;
    Eigen::MatrixXd *xData;
    Eigen::VectorXd *yData;
    set<int> constraintSet;
    
    double th;    
    state *nextRandomStateDown();
    state *nextRandomStateUp();
    state *nextUnexploredStateDown();
    state *nextUnexploredStateUp();
            
  public:
    static int totalIter;
    static unordered_map<std::string, state*> stateMap;
    static int maxInls;
    static Eigen::VectorXd bestTheta;
    static pType problemType;
    static Eigen::MatrixXd *AA;
    static Eigen::MatrixXd *bb;
    static Eigen::MatrixXd *cc;
    static Eigen::VectorXd *dd;
    
    double minimax;
    double state_th;
    //vec theta;
    Eigen::VectorXd xn;    
    stateDirection direction;
    set<int> basics;
    set<int> violationSet;
    set<int> supportSet;    
    
    int inls;
    int ptsInls;
    state();
    state(Eigen::MatrixXd *xData, Eigen::VectorXd *yData, const set<int>&vSet, const set<int> &sSet, double th_, stateDirection prevDirection, const Eigen::VectorXd &theta0);    
    //~state();    
    
    int getTotalIter();
    
    state *nextRandomState();
    state *nextUnexploredState();

    bool terminated();
    bool feasible();
    int computeReward();
    int nInliers();


  };

   
}


#endif


//    state(mat *x_, vec *y_, const set<int>&vSet, const set<int> &sSet, double th_, stateDirection prevDirection, const vec &theta0);

