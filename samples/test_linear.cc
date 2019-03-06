#include <iostream>
#include <string>
#include <cstring>

#include <glog/logging.h>
#include <gflags/gflags.h>

#include <cstdlib>
#include <ctime>
#include <set>
#include <Eigen/Core>

#include "problems/LinearEstimator.hh"
#include "common.hh"
#include "usac.hh"
#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "State.hh"
#include "tree.hh"
#include "mcts.hh"
#include "LinearFit.hh"


using namespace std;
using namespace usac;
using namespace arma;
using namespace LinearTools;
using namespace MCTS;
using namespace ToolBox;

DEFINE_bool(use_prosac, false, "");
DEFINE_int32(N,100, "Data size");
DEFINE_int32(d,8, "Data dimensions");
DEFINE_double(th, 0.06, "Inlier Threshold");
DEFINE_double(p, 0.2, "Outlier Percentage");

int main(int argc, char **argv){
  
  std::string FLAGS_config;
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  LinearEstimator sac;  
  ConfigFileReader cfr;
  LinearParams linear_params;
  USACParams usac_params;

  int nr;
  int nOutlier;
  int d;
  Eigen::MatrixXd xData;
  Eigen::VectorXd yData;
  Eigen::VectorXd thetaD;
  mat x;
  vec y;
  vec theta;
  double minimax;
  const double insig = 0.01;
  const double osig = 0.1;
  
  double inlier_th;
  
  nr = FLAGS_N;
  
  double outlierP = FLAGS_p;
  cout <<outlierP<<"\n" ;
  nOutlier = round(nr*outlierP);
  d = FLAGS_d;
  
  cfr.readFile(argv[1]);
  CHECK(readParamsFromConfigFile(cfr, &usac_params,&inlier_th));
  usac_params.min_sample_size = d;
      
  inlier_th = FLAGS_th;

  //Gen random data  
  srand(time(0));
  genRandomData(nr, d, insig, osig, outlierP, xData, yData, thetaD);
  
  
  LinearFit linearFitter(&xData, &yData, inlier_th);
  linearFitter.configFile = argv[1];
  linearFitter.sampleFile = "samples.txt";
    
  
  linearFitter.Fit("RANSAC");  
  //linearFitter.writeSamplesToFile("samples.txt");  
  cout << "Inls = " << linearFitter.inls << " Run Time = " << linearFitter.runTime << "\n";;

  linearFitter.Fit("MCTS");
  cout << "Inls = " << linearFitter.inls << " Run Time = " << linearFitter.runTime << "\n";;
  
  
  linearFitter.Fit("ASTAR");    
  cout << "Inls = " << linearFitter.inls << " Run Time = " << linearFitter.runTime << "\n";;
  
  
   
  
  //linearFitter.Fit("LORANSAC");
  //cout << "Inls = " << linearFitter.inls << " Run Time = " << linearFitter.runTime << "\n";;
  
  
  

  // //SAC Solution
  // sac.allSamples.clear();
  // sac.init_params(usac_params);
  // sac.init_data(nr, inlier_th);
  
  // sac.init_problem(linear_params, &xData, &yData);
  // sac.solve();
  // print_usac_results(sac.results_);
  // sac.writeSamplesToFile("samples.txt");

  
  
  // sac.sampling_method = SAMP_FILE;
  // sac.init_data(nr, inlier_th);
  // sac.init_problem(linear_params, &xData, &yData);
  // sac.readSamplesFromFile("samples.txt"); 
  // sac.solve();
  // print_usac_results(sac.results_);
  
  
  return 0;
}

	   


//

// set<int> basics;
// Eigen::VectorXd theta0f = Eigen::VectorXd::Random(d);
// Eigen::VectorXd thetaf = Chebyshev::Descend(xData, yData, inlier_th, theta0f, basics, minimax);
// int inls = findInliers(xData, yData, thetaf, inlier_th);
  

// set<int>vSet;  
// set<int>sSet;

// vSet.clear();
// for (int i=0; i<nr; i++)
//   sSet.insert(i);

// state::totalIter = 0;
// state::maxInls = 0;
// state::stateMap.clear();

  
  // clock_t tick = clock();
  // state initState(&xData, &yData, vSet, sSet, inlier_th, MCTS::DOWN, thetaf);
  // //Create a tree with the init state;
  // Tree *searchTree = new Tree(&initState);
  // Node *currentNode = searchTree->root;  
  // while (!currentNode->nodeState->terminated()){
  //   state *currentState = currentNode->nodeState;    
  //   cout << currentState->minimax <<endl;// << "  "<<currentNode<<endl;
  //   int maxChildren = currentState->direction==DOWN?
  //     currentState->basics.size():currentState->violationSet.size();    
  //   int iter = maxChildren;    
  //   Node *prevNode = currentNode;
  //   currentNode = uctSearch(searchTree, currentNode, iter);
    
  //   if (prevNode == currentNode){
  //     cout << "---------=======-------------"<<endl;
  //     break;
  //   }
  // }


  // double mcts_time = (clock()-tick)*1.0/CLOCKS_PER_SEC;
  
  //inls = findInliers(xData, yData, currentNode->nodeState->xn, inlier_th);
  //inls = state::maxInls;
  
  
  //cout << "MCTS Inls : "<<inls <<endl;
  //  cout << "MCTS Time: " <<mcts_time << endl;
  //  cout << "MCTS Iter: " <<state::totalIter << endl;

