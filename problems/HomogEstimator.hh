#ifndef HOMOESTIMATOR_HH
#define HOMOESTIMATOR_HH

#include "usac.hh"
#include <Eigen/Core>
#include "quasiconvexFunctions.hh"

using namespace Eigen;

namespace usac {

struct HomogEstimator: public USAC
{
  ~HomogEstimator();

  bool init_problem(int nr, double* pointData);
  void cleanup_problem();
  double compute_error(int model_idx, int sample_idx);
  int gen_min_sample_models(const vector<int> &sample);
  bool gen_refined_model(const vector<int>& sample, const vector<double> *weights = NULL);
  bool validate_sample(const vector<int> &sample);
  void store_model(int modelIndex, int numInliers);

  vector<double> final_model_params_;

  /*Matrices for computing L1 residual*/
  MatrixXd *AA, *bb, *cc;
  VectorXd *dd;
  /*Matrices for Linear Constraints*/
  MatrixXd *xData;
  VectorXd *yData;
    
 private:
  double* input_points_denorm_ = NULL;
  double* input_points_ = NULL;
  double* data_matrix_ = NULL;
  double  m_T1_[9], m_T2_[9], m_T2inv_[9];
  vector<double*> models_;
  vector<double*> models_denorm_;
  double inv_model[9];
};

} // namespace usac

#endif
