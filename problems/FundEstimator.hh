#ifndef FUNDESTIMATOR_HH
#define FUNDESTIMATOR_HH

#include "usac.hh"

namespace usac {

  enum MatrixDecomposition      {DECOMP_QR, DECOMP_LU};
  struct FundParams
  {
    MatrixDecomposition decomp_alg = DECOMP_QR;
    double h_degen_th = .0;
    int max_upgrade_samples = 500;
  };

  class FundEstimator: public USAC
  {
  public:
    bool init_problem(const FundParams& cfg, int nr, double* pointData);
    // ------------------------------------------------------------------------
    // storage for the final result
    vector<double> final_model_params_;
    vector<double> degen_final_model_params_;

  public:
    FundEstimator()
    {
      input_points_ = NULL;
      data_matrix_  = NULL;
      degen_data_matrix_  = NULL;
      models_.clear();
      models_denorm_.clear();
    };
    ~FundEstimator()
    {
      if (input_points_) { delete[] input_points_; input_points_ = NULL; }
      if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
      if (degen_data_matrix_) { delete[] degen_data_matrix_; degen_data_matrix_ = NULL; }
      for (size_t i = 0; i < models_.size(); ++i)
	{
	  if (models_[i]) { delete[] models_[i]; }
	}
      models_.clear();
      for (size_t i = 0; i < models_denorm_.size(); ++i)
	{
	  if (models_denorm_[i]) { delete[] models_denorm_[i]; }
	}
      models_denorm_.clear();
    };

  public:
    void cleanup_problem();
    int gen_min_sample_models(const vector<int> &sample);
    bool gen_refined_model(const vector<int>& sample, const vector<double> * weights = NULL);
    bool validate_model(int modelIndex, const vector<int> &sample);
    double compute_error(int model_idx, int sample_idx);
    void test_soln_degeneracy(bool* degenerateModel, bool* upgradeModel);
    virtual int upgrade_degenerate_model(vector<double> *errs);
    void find_weights(int model_idx, const vector<int>& inliers, vector<double>*weights);
    void find_weights(int modelIndex, const vector<int> &inliers,
		      int numInliers, double* weights);
    void store_model(int modelIndex, int numInliers);

  private:
    double* input_points_denorm_;
    double degen_homog_threshold_ = .0;
    int degen_max_upgrade_samples_ = 0;
    MatrixDecomposition matrix_decomposition_method;
    double* input_points_;
    double* data_matrix_;
    double* degen_data_matrix_ = NULL;
    vector<int> degen_outlier_flags_;
    double  m_T1_[9], m_T2_[9], m_T2_trans_[9];
    vector<double*> models_;
    vector<double*> models_denorm_;
  };

} //namespace usac

#endif
