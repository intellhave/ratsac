#ifndef LINEARESTIMATOR_HH
#define LINEARESTIMATOR_HH

#include "usac.hh"
//#include <armadillo>
#include <Eigen/Core>
#include <Eigen/Dense>

//using namespace arma;
using namespace std;
using namespace Eigen;
struct LinearParams
{
  
};


namespace usac {
  class LinearEstimator: public USAC{
  private:
    int npts;
    int d;
    Eigen::MatrixXd* xData;
    Eigen::VectorXd* yData;
    vector<double*> models_;
    vector<double*> models_denorm_;
       
  public:
    bool init_problem(const LinearParams& cfg, Eigen::MatrixXd *x_, Eigen::VectorXd* y_);
    Eigen::VectorXd final_model_params_;
    //vec degen_final_model_params_;
    
  public:
    LinearEstimator(){
    }

  public:
    void cleanup_problem();
    int gen_min_sample_models(const vector<int>&sample);
    double compute_error(int model_idx, int sample_idx);
    void store_model(int modelIndex, int numInliers);
    
    bool gen_refined_model(const vector<int>&sample, const vector<double>*weights = NULL);
    bool validate_model(int modelIndex, const vector<int>&sample);
    void find_weights(int model_idx, const vector<int>&inliers, vector<double>* weights);
    
		      
  };
  
}




#endif
