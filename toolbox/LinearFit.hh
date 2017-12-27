#ifndef LINEARFIT_HH
#define LINEARFIT_HH


#include <cstdlib>
#include <Eigen/Core>
#include <vector>
#include <set>
#include <limits>
#include "common.hh"
#include <ctime>

#include "LinearEstimator.hh"
#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "State.hh"
#include "tree.hh"
#include "mcts.hh"
#include "usac.hh"
#include "Astar.hh"
#include "extra.hh"
#include "config.hh"
#include "LinearProgramming.hh"

using namespace std;
using namespace LinearTools;
using namespace MCTS;
using namespace usac;

namespace ToolBox{

class LinearFit{
private:
  double inlier_th;                //Inlier Threshold
  Eigen::MatrixXd *xData;  
  Eigen::VectorXd *yData;
  LinearEstimator sac;
  
  clock_t startTime;
  clock_t stopTime;
  int N;
  int d;
  set<int> basics;
  void ransacFit();
  void MCTSFit();
  void L1Fit();
  void LInfFit();
  
  
  ConfigFileReader cfr;

  Config *cfg;  
  
  RandomSamplingMethod sampling_method;
  VerifMethod verif_method;
  bool need_local_optim;
  bool mle_sac;

  
  
public:
  LinearFit(Eigen::MatrixXd *x_, Eigen::VectorXd *y_, Config *cfg_); //Constructor
  //~LinearFit();
  Eigen::VectorXd bestTheta;  //Fitting Results
  int inls; //Number of Inliers;
  std::vector<int> consensusSet;
  double runTime;
  const char *configFile;
  std::string sampleFile;

  std::vector<int> *sorted;
  
  void Fit(const std::string& method);
  
  

  
  
  
  
  
};


}



#endif

