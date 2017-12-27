#ifndef QUASICONVEX_HH
#define QUASICONVEX_HH

#include <cstdlib>
#include <Eigen/Core>
#include <vector>
#include <set>
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
#include "common.hh"
#include "problems/HomogEstimator.hh"
#include "problems/AffineEstimator.hh"
#include "config.hh"
#include "quasiconvexFunctions.hh"
#include "MSearch.hh"
#include "extra.hh"
#include "LinearProgramming.hh"

using namespace std;
using namespace LinearTools;
using namespace MCTS;
using namespace usac;

namespace ToolBox{
  

  
  class QuasiconvexFit{
  private:
    /*Input Data*/
    double inlier_th;
    Eigen::MatrixXd u1;
    Eigen::MatrixXd u2;
    std::vector<int> *sorted; /*Sorted indexes based on matching scores*/

    /*Data used for Chebyshev fit (L1, Linf)*/
    Eigen::MatrixXd xData;
    Eigen::VectorXd yData;
    /*Quasiconvex matrices*/
    Eigen::MatrixXd A, b, c;
    Eigen::VectorXd d;

    clock_t startTime;
    clock_t stopTime;

    int nPoints;
    int N;
    int dim;
    set<int> basics;
    /*General Configurations*/
    Config *cfg;    
    /*For RANSAC variants;*/
    ConfigFileReader cfr;
    RandomSamplingMethod sampling_method;
    VerifMethod verif_method;
    bool need_local_optim;
    bool mle_sac;
    
    prob problem;
    void ransacFit();
    void affineRansacFit();
    void homoRansacFit();    
    void MCTSFit();
    void genMatrixHomography(); /*From input points, genearte A, b, c,d for quasiconvex residual*/
    void genMatrixAffine();
    void genMatrixQuasiconvex();
    
    void LInfFit();
    void L1Fit();


  public:
    QuasiconvexFit(const Eigen::MatrixXd &u1_, const Eigen::MatrixXd &u2_, std::vector<int> *sorted_, Config *cfg_);
    Eigen::VectorXd bestTheta;
    int inls;

    double runTime;
    const char *configFile;
    std::string sampleFile;

    void Fit(const std::string& method);


  };


}






#endif
