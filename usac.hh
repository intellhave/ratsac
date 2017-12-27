#ifndef USAC_HH
#define USAC_HH

#include <glog/logging.h>
#include <time.h>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <list>
#include <functional>
#include <memory>
#include <fstream>
#include <cstring>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <math.h>

namespace usac {

  using std::vector;
  using std::list;
  using std::shared_ptr;
  using std::unique_ptr;

  template<typename T>
  static inline void resize_fill_vec(vector<T> &vec, int n, T v=0) {
    vec.resize(n);
    if (!v)
      std::fill(vec.begin(), vec.end(), v);
  }

  enum RandomSamplingMethod {SAMP_UNIFORM, SAMP_PROSAC, SAMP_FILE};
  enum VerifMethod {VERIF_STANDARD, VERIF_SPRT};

  extern const double USAC_CONF_TH;
  extern const int USAC_MAX_HYPS;
  extern const int USAC_PROSAC_MAX_SAMPLES;
  extern const int USAC_MIN_SAMPLES;
  extern const double USAC_PROSAC_BETA;
  extern const double USAC_PROSAC_NON_RAND_CONF;
  extern const int USAC_MIN_STOP_LEN;

  extern const double USAC_SPRT_tM;
  extern const double USAC_SPRT_mS;
  extern const double USAC_SPRT_DELTA;
  extern const double USAC_SPRT_EPS;

  extern const int USAC_LOSAC_INNER_SAMPLE_SIZE;
  extern const double USAC_LOSAC_INNER_SAMPLE_SIZE_RATIO;
  extern const int USAC_LOSAC_INNER_RANSAC_REPS;
  extern const double USAC_LOSAC_TH_MULTIPLIER;
  extern const int USAC_LOSAC_ITER_STEPS;

  struct USACParams
  {
    double conf_th = USAC_CONF_TH;
    int min_sample_size = 0;
    int min_inls = 0;
    int max_hypotheses = USAC_MAX_HYPS;
    int min_hypotheses = USAC_MIN_SAMPLES;
    int max_solns_per_sample = 1;
    bool need_prevalidate_sample = false;
    bool need_prevalidate_model = false;
    bool need_test_degeneracy = false;
    bool need_local_optim = false;
    bool mle_sac = false;
    
    RandomSamplingMethod sampling_method = SAMP_UNIFORM;
    VerifMethod verif_method = VERIF_STANDARD;

    int  prosac_max_samples = USAC_PROSAC_MAX_SAMPLES;
    double prosac_beta = USAC_PROSAC_BETA;
    double prosac_non_rand_conf = USAC_PROSAC_NON_RAND_CONF;
    int  prosac_min_stop_len = USAC_MIN_STOP_LEN;

    double sprt_tM = USAC_SPRT_tM;
    double sprt_mS = USAC_SPRT_mS;
    double sprt_delta = USAC_SPRT_DELTA;
    double sprt_epsilon = USAC_SPRT_EPS;

    int losac_inner_sample_size = USAC_LOSAC_INNER_SAMPLE_SIZE;
    int losac_inner_ransac_reps = USAC_LOSAC_INNER_RANSAC_REPS;
    double losac_threshold_multiplier = USAC_LOSAC_TH_MULTIPLIER;
    int losac_num_iterative_steps = USAC_LOSAC_ITER_STEPS;
  };

  struct USAC;

  struct Prosac
  {
    Prosac() {}
    Prosac(int max_samples_, double beta_, double non_rand_conf_, int min_stop_len_): max_samples(max_samples_), beta(beta_), non_rand_conf(non_rand_conf_), min_stop_len(min_stop_len_) {}

    void init_data(int data_size, int min_sample_size, int max_hypotheses, const vector<int> &sorted_indices);
    void gen_min_sample(USAC *u, vector<int> *sample);
    int update_stopping(USAC *u, int hyp_nr);

    int  max_samples = USAC_PROSAC_MAX_SAMPLES;
    double beta = USAC_PROSAC_BETA;
    double non_rand_conf = USAC_PROSAC_NON_RAND_CONF;
    int  min_stop_len = USAC_MIN_STOP_LEN;

    int subset_size_;
    int largest_size_;
    int stop_len_;
    vector<int> growth_function_;
    vector<int> non_random_inliers_;
    vector<int> maximality_samples_;
    vector<int>  sorted_indices_;
  };

  struct Sprt
  {
    Sprt() {}
    Sprt(double tM_, double mS_, double delta_, double epsilon_): tM(tM_), mS(mS_), delta(delta_), epsilon(epsilon_) {}

    void init_data();

    double design_test(double delta, double eps);
    void design_cur_test();
    void add_test_history(int hyp_nr);
    bool evaluate_model(USAC *u, int model_idx, vector<double> *errs, int* inlier_nr, int* test_nr);
    bool verify_and_store(USAC *u, int model_idx, int inlier_nr, int test_nr, bool good);      
    void update_stopping(USAC *u, int inl_nr, int *stopping_nr);

    double tM = USAC_SPRT_tM;
    double mS = USAC_SPRT_mS;
    double delta = USAC_SPRT_DELTA;
    double epsilon = USAC_SPRT_EPS;

    double cur_delta_;
    double cur_epsilon_;
    double decision_threshold_;
    bool need_update_stopping_;
    int last_wald_history_update_;

    struct TestHistorySPRT
    {
      double epsilon_;
      double delta_;
      double A_;
      int k_;
    };
    list<TestHistorySPRT> wald_test_history_;
  };

  struct Losac
  {
    Losac() {}
    Losac(int inner_sample_size_, int inner_ransac_reps_, double threshold_multiplier_, int num_iterative_steps_): inner_sample_size(inner_sample_size_), inner_ransac_reps(inner_ransac_reps_), threshold_multiplier(threshold_multiplier_), num_iterative_steps(num_iterative_steps_) {}

    void init_data(int data_size) {
      resize_fill_vec(prev_best_inliers_, data_size, false);
      num_prev_best_inliers_ = 0;
    }

    int locally_optimize_soln(USAC *u, int best_inliers);
    void optimize_and_update(USAC *u);

    int inner_sample_size = USAC_LOSAC_INNER_SAMPLE_SIZE;
    int inner_ransac_reps = USAC_LOSAC_INNER_RANSAC_REPS;
    double threshold_multiplier = USAC_LOSAC_TH_MULTIPLIER;
    int num_iterative_steps = USAC_LOSAC_ITER_STEPS;

    int num_prev_best_inliers_;
    vector<bool> prev_best_inliers_;
  };

  struct USACResults {
    void reset() {
      hyp_nr_ = 0;
      model_nr_ = 0;
      rejected_sample_nr_ = 0;
      rejected_model_nr_ = 0;
      best_inlier_nr_ = 0;
      degen_inlier_nr_ = 0;
      total_points_verified_ = 0;      
      num_local_optimizations_ = 0;
      min_mle_penalty_ = std::numeric_limits<double>::max();      
      inlier_flags_.clear();
      best_sample_.clear();
      degen_inlier_flags_.clear();
      degen_sample_.clear();
    }
    void resize(int nr, int min_sample_size, bool need_test_degeneracy) {
      data_size = nr;
      inlier_flags_.resize(nr, 0);
      best_sample_.resize(min_sample_size, 0);
      if (need_test_degeneracy) {
	degen_inlier_flags_.resize(nr, 0);
	degen_sample_.resize(min_sample_size, 0);
      }
    }
    int data_size = 0;
    int hyp_nr_ = 0;
    int model_nr_ = 0;
    int rejected_sample_nr_ = 0;
    int rejected_model_nr_ = 0;
    int total_points_verified_ = 0;
    int num_local_optimizations_ = 0;
    int best_inlier_nr_ = 0;
    double min_mle_penalty_ = std::numeric_limits<double>::max();
    vector<bool> inlier_flags_;
    vector<int> best_sample_;

    int degen_inlier_nr_ = 0;
    vector<bool> degen_inlier_flags_;
    vector<int> degen_sample_;
  };

  enum IterState {ITER_END, ITER_CONTINUE, ITER_NOT_ENOUGH_DATA};
  
  int find_inliers(const vector<double>& errs, double threshold, vector<int>* inliers);
  void gen_uniform_rand_sample(int dataSize, int sampleSize, vector<int>* sample);
  void get_sample_from_list(const vector<unsigned int>& sampleList, unsigned int &currIdx, int sampleSize, vector<int>* sample);
  
  struct USAC
  {
    void init_params(const USACParams &p); /* run once until params changes */
    void init_data(int sample_nr, double inlier_th, const vector<int> *prosac_sorted_indices=NULL); /* run before new data */
    void writeSamplesToFile(const std::string& fileName);
    void readSamplesFromFile(const std::string& fileName);
    unsigned int currentSampleIdx;
    
    bool solve();

    IterState iterate(int it_nr);

    USACResults results_;

    virtual int gen_min_sample_models(const vector<int> &sample) = 0;
    virtual bool gen_refined_model(const vector<int>& sample, const vector<double> *weights=NULL) = 0;
    virtual bool validate_sample(const vector<int>&) {return true;}
    virtual bool validate_model(int, const vector<int> &) {return true;}
    bool evaluate_model_ransac(int model_idx, vector<double> *errs, int *inlier_nr, int *test_nr);
    virtual double compute_error(int model_index, int sample_idx) {return 0;}
    virtual void test_soln_degeneracy(bool* degenerateModel, bool* upgradeModel) {*degenerateModel = false; *upgradeModel = false;};
    virtual int upgrade_degenerate_model(vector<double> *errs) {return 0;};
    virtual void find_weights(int model_idx, const vector<int>& inliers,
			      vector<double> *weights) {
      resize_fill_vec(*weights, inliers.size(), 1.0);
    };    
    virtual void store_model(int model_idx, int numInliers) {};

    bool verify_and_store(int model_idx, int inlier_nr, int test_nr, bool good);
    bool verify_and_store_mlesac(int model_idx, double penalty, int inlier_nr);
				 
    /*For MLESAC*/
    double penaltyMLE(vector<double> *errs, double dSigma);
    
    double conf_th = USAC_CONF_TH;
    int min_sample_size = 0;
    int min_inls = 0;
    int max_hypotheses = USAC_MAX_HYPS;
    int min_hypotheses = USAC_MIN_SAMPLES;
    int max_solns_per_sample = 1;
    int nEMIter = 3;
    
    bool need_prevalidate_sample = false;
    bool need_prevalidate_model = false;
    bool need_test_degeneracy = false;
    bool need_local_optim = false;
    bool mle_sac = false;
    
    RandomSamplingMethod sampling_method = SAMP_UNIFORM;
    VerifMethod verif_method = VERIF_SPRT;

    unique_ptr<Prosac> prosac = NULL;
    unique_ptr<Sprt> sprt = NULL;
    unique_ptr<Losac> losac = NULL;

    std::function<bool(int, vector<double>*, int*, int*)> evaluate_model;

    int data_size_;
    double inlier_th_;
    vector<int>  min_sample_;
    vector<double> tmp_errs_;
    vector<double>  best_errs_;

    int adaptive_stopping_nr_;

    int update_stopping(int inl_nr, int total_pts, int sample_size);
    void store_soln(int model_idx, int inl_nr);

    static std::vector<unsigned int>allSamples;

  
  };

} // namespace usac
#endif
