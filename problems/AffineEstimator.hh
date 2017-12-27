#ifndef AFFINE_HH
#define AFFINE_HH

#include "usac.hh"
#include "quasiconvexFunctions.hh"
#include <Eigen/Core>
#include <Eigen/Dense>

struct AffineParams{
  
};

namespace usac{
  
  class AffineEstimator: public USAC{
  private:
    int npts;
    Eigen::MatrixXd *u1;
    Eigen::MatrixXd *u2;
    Eigen::MatrixXd *AA, *bb, *cc;
    Eigen::VectorXd *dd;                    
    vector<double*> models_;

  public:
    bool init_problem(Eigen::MatrixXd *u1, Eigen::MatrixXd *u2);
    bool init_problem(Eigen::MatrixXd *A, Eigen::MatrixXd *b, Eigen::MatrixXd *c, Eigen::VectorXd *d);
        
    
    Eigen::VectorXd final_model_params_;

    AffineEstimator(){
    }

    void cleanup_problem();
    int gen_min_sample_models(const vector<int>& sample);
    double compute_error(int model_idx, int sample_idx);
    void store_model(int modelIndex, int numInliers);
    bool gen_refined_model(const vector<int> &sample, const vector<double>* weights = NULL);
    bool validate_model(int modelIndex, const vector<int>& sample);
    void find_weights(int model_idx, const vector<int>& inliers, vector<double>* weights);
    
  };
  
}


#endif
