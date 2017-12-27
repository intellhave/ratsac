#include <math.h>
#include "usac.hh"

namespace usac {

  using namespace std::placeholders;
  static const double EPS = std::numeric_limits<double>::epsilon();

  const double USAC_CONF_TH = 0.95;
  const int USAC_MAX_HYPS = 1000;
  const int USAC_MIN_SAMPLES = 100;
  const int USAC_PROSAC_MAX_SAMPLES = 2000;
  const double USAC_PROSAC_BETA = 0.05;
  const double USAC_PROSAC_NON_RAND_CONF = 0.95;
  const int USAC_MIN_STOP_LEN = 20;

  const double USAC_SPRT_tM = 200.0;
  const double USAC_SPRT_mS = 2.38;
  const double USAC_SPRT_DELTA = 0.05;
  const double USAC_SPRT_EPS = 0.2;

  const int USAC_LOSAC_INNER_SAMPLE_SIZE = 10;
  const double USAC_LOSAC_INNER_SAMPLE_SIZE_RATIO = 0.9;
  const int USAC_LOSAC_INNER_RANSAC_REPS = 10;
  const double USAC_LOSAC_TH_MULTIPLIER = 2.0;
  const int USAC_LOSAC_ITER_STEPS = 4;

  int find_inliers(const vector<double>& errs, double threshold, vector<int>* inliers)
  {
    inliers->clear();
    for (int i = 0; i < (int) errs.size(); i++)
      if ((errs[i] < threshold) && (errs[i]>0))
	inliers->push_back(i);
    return inliers->size();
  }

  void gen_uniform_rand_sample(int dataSize, int sampleSize, vector<int>* sample)
  {
    srand(time(NULL));
    sample->resize(sampleSize);
    int count = 0;
    int index;
    auto pos = sample->begin();
    do {
      index = rand() % dataSize;
      if (find(sample->begin(), pos, index) == pos) {
	(*sample)[count] = index;
	++count;
	++pos;
      }
    } while (count < sampleSize);
  }

  void get_sample_from_list(const vector<unsigned int>& sampleList, unsigned int &currIdx, int sampleSize, vector<int>* sample){
    assert(currIdx <sampleList.size());
    (*sample).clear();
    for (int i=0; i<sampleSize; i++){
      (*sample).push_back(sampleList[currIdx]);
      currIdx++;      
    }       
  }
  
  void Prosac::init_data(int data_size, int min_sample_size, int max_hypotheses, const vector<int> &sorted_indices)
  {
    CHECK(sorted_indices.size());
    sorted_indices_ = sorted_indices;
    growth_function_.clear();
    growth_function_.resize(data_size);
    
    double T_n;
    int T_n_p = 1;
    // compute initial value for T_n
    T_n = max_samples;
    for (int i = 0; i < min_sample_size; ++i)
      T_n *= (min_sample_size-i) * 1.0/(data_size-i);

    // compute values using recurrent relation
    for (int i = 0; i < data_size; ++i) {
      if (i+1 <= min_sample_size) {
	growth_function_[i] = T_n_p;
	continue;
      }
      double temp = (double)(i+1)*T_n / (i+1-min_sample_size);
      growth_function_[i] = T_n_p + (int)ceil(temp - T_n);
      T_n = temp;
      T_n_p = growth_function_[i];
    }

    // ------------------------------------------------------------------------
    // initialize the data structures that determine stopping

    // non-randomness constraint
    // i-th entry - inlier counts for termination up to i-th point (term length = i+1)
    non_random_inliers_.clear();
    non_random_inliers_.resize(data_size, 0);
    double pn_i = 1.0;    // prob(i inliers) with subset size n
    for (int n = min_sample_size+1; n <= data_size; ++n) {
      if (n-1 > 1000) {
	non_random_inliers_[n-1] = non_random_inliers_[n-2];
	continue;
      }

      vector<double> pn_i_vec(data_size, 0);
      // initial value for i = m+1 inliers
      pn_i_vec[min_sample_size] = (beta)*std::pow((double)1-beta, (double)n-min_sample_size-1)*(n-min_sample_size);
      pn_i = pn_i_vec[min_sample_size];
      for (int i = min_sample_size+2; i <= n; ++i) {
	// use recurrent relation to fill in remaining values
	if (i == n) {
	  pn_i_vec[n-1] = std::pow((double)beta, (double)n-min_sample_size);
	  break;
	}
	pn_i_vec[i-1] = pn_i * ((beta)/(1-beta)) * ((double)(n-i)/(i-min_sample_size+1));
	pn_i = pn_i_vec[i-1];
      }
      // find minimum number of inliers satisfying the non-randomness constraint
      double acc = 0.0;
      int i_min = 0;
      for (int i = n; i >= min_sample_size+1; --i) {
	acc += pn_i_vec[i-1];
	if (acc < 1-non_rand_conf)
	  i_min = i;
	else
	  break;
      }
      non_random_inliers_[n-1] = i_min;
    }

    resize_fill_vec(maximality_samples_, data_size, max_hypotheses);
    largest_size_ = min_sample_size;
    subset_size_ = min_sample_size;
    stop_len_    = data_size;
  }

  void Prosac::gen_min_sample(USAC *u, vector<int>* sample)
  {
    int hypCount = u->results_.hyp_nr_;
    if (hypCount > max_samples) {
      gen_uniform_rand_sample(u->data_size_, u->min_sample_size, sample);
      return;
    }

    // if current stopping length is less than size of current pool, use only points up to the stopping length
    if (subset_size_ > stop_len_) {
      gen_uniform_rand_sample(stop_len_, u->min_sample_size, sample);
      return;
    }
    
    // increment the size of the sampling pool if required
    if (hypCount > growth_function_[subset_size_-1]) {
      ++subset_size_;
      if (subset_size_ > u->data_size_)
	subset_size_ = u->data_size_;

      if (largest_size_ < subset_size_)
	largest_size_ = subset_size_;
    }

    // generate PROSAC sample
    gen_uniform_rand_sample(subset_size_-1, u->min_sample_size-1, sample);
    (*sample)[u->min_sample_size-1] = subset_size_-1;
    for (size_t i = 0; i < sample->size(); ++i)
      (*sample)[i] = sorted_indices_[(*sample)[i]];
  }

  int Prosac::update_stopping(USAC *u, int hyp_nr)
  {
    if (hyp_nr > max_samples)
      return u->update_stopping(u->results_.best_inlier_nr_, u->data_size_, u->min_sample_size);

    int max_samples = maximality_samples_[stop_len_-1];

    int inlier_count = 0;
    for (int i = 0; i < min_stop_len; ++i)
      inlier_count += u->results_.inlier_flags_[sorted_indices_[i]];

    // after this initial subset, try to update the stopping length if possible
    for (int i = min_stop_len; i < u->data_size_; ++i) {
      inlier_count += u->results_.inlier_flags_[sorted_indices_[i]];

      if (non_random_inliers_[i] < inlier_count) {
	non_random_inliers_[i] = inlier_count;   // update the best inliers for the the subset [0...i]

	// update the number of samples based on this inlier count
	if ( (i == u->data_size_ - 1) ||
	     (u->results_.inlier_flags_[sorted_indices_[i]] && !u->results_.inlier_flags_[sorted_indices_[i+1]]))
	  {
	    int new_samples = u->update_stopping(inlier_count, i+1, u->min_sample_size);
	    if (i+1 < largest_size_)
	      new_samples += hyp_nr - growth_function_[i];

	    if (new_samples < maximality_samples_[i]) {
	      // if number of samples can be lowered, store values and update stopping length
	      maximality_samples_[i] = new_samples;
	      if ( (new_samples < max_samples) || ( (new_samples == max_samples) && (i+1 >= stop_len_) ) )
		{
		  stop_len_ = i+1;
		  max_samples = new_samples;
		}
	    }
	  }
      }
    }
    return max_samples;
  }

  void Sprt::init_data()
  {
    cur_delta_ = delta;
    cur_epsilon_ = epsilon;
    last_wald_history_update_ = 0;
    wald_test_history_.clear();
    design_cur_test();
    need_update_stopping_ = true;
  }

  double Sprt::design_test(double cur_delta, double cur_eps)
  {
    double An_1, An, C, K;

    C = (1 - cur_delta)*log( (1 - cur_delta)/(1-cur_eps) )
      + cur_delta*(log( cur_delta/cur_eps));
    K = (tM*C)/mS + 1;
    An_1 = K;

    // compute A using a recursive relation
    // A* = lim(n->inf)(An), the series typically converges within 4 iterations
    for (int i = 0; i < 10; ++i) {
      An = K + log(An_1);
      if (An - An_1 < 1.5e-8)
	break;
      An_1 = An;
    }
    return An;
  }

  void Sprt::design_cur_test()
  {
    decision_threshold_ = design_test(cur_delta_, cur_epsilon_);
  }

  void Sprt::add_test_history(int hyp_nr)
  {
    TestHistorySPRT new_th;
    new_th.epsilon_ = cur_epsilon_;
    new_th.delta_ = cur_delta_;
    new_th.A_ = decision_threshold_;
    new_th.k_ = hyp_nr - last_wald_history_update_;
    last_wald_history_update_ = hyp_nr;

    wald_test_history_.push_front(new_th);
  }

  bool Sprt::evaluate_model(USAC *u, int model_idx, vector<double> *errs, int* inlier_nr, int* test_nr)
  {
    bool good_flag = true;
    double lambdaj, lambdaj_1 = 1.0;
    int inl_nr = 0;
    int i;

    for (i = 0; i < (int) errs->size(); ++i) {
      double err = u->compute_error(model_idx, i);
      (*errs)[i] = err;
      if (err < u->inlier_th_) {
	inl_nr++;
	lambdaj = lambdaj_1 * (cur_delta_/cur_epsilon_);
      } else {
	lambdaj = lambdaj_1 * ((1 - cur_delta_)/(1 - cur_epsilon_));
      }

      if (lambdaj > decision_threshold_) {
	good_flag = false;
	break;
      } else {
	lambdaj_1 = lambdaj;
      }
    }

    if (inlier_nr)
      *inlier_nr = inl_nr;
    if (test_nr)
      *test_nr = i + 1;
    return good_flag;
  }

  bool Sprt::verify_and_store(USAC *u, int model_idx, int inlier_nr, int test_nr, bool good)
  {
    if (!good) {
      double delta_new = inlier_nr * 1.0 / test_nr;
      if (delta_new > 0 && abs(cur_delta_ - delta_new) / cur_delta_ > 0.1) {
	// update parameters
	add_test_history(u->results_.hyp_nr_);
	design_cur_test();
	cur_delta_ = delta_new;
      }
    } else if (inlier_nr > u->results_.best_inlier_nr_ && inlier_nr > u->min_inls) {
      need_update_stopping_ = true;
      u->results_.best_inlier_nr_ = inlier_nr;

      add_test_history(u->results_.hyp_nr_);
      cur_epsilon_ = u->results_.best_inlier_nr_ * 1.0 / u->data_size_;
      design_cur_test();
      u->store_soln(model_idx, u->results_.best_inlier_nr_);
      return true;
    }
    return false;
  }

  static double compute_exp_sprt(double newEpsilon, double epsilon, double delta)
  {
    double al, be, x0, x1, v0, v1, h;

    al = log(delta/epsilon);
    be = log( (1-delta)/(1-epsilon) );

    x0 = log( 1/(1-newEpsilon) )/be;
    v0 = newEpsilon * exp(x0 *al);
    x1 = log( (1-2*v0) / (1-newEpsilon) )/be;
    v1 = newEpsilon * exp(x1 * al) + (1-newEpsilon) * exp(x1 * be);
    h = x0 - (x0 - x1)/(1+v0 - v1)*v0;
    return h;
  }

  void Sprt::update_stopping(USAC *u, int inl_nr, int *stopping_nr)
  {
    if (!(u->results_.hyp_nr_ >= *stopping_nr && need_update_stopping_))
      return;

    double n_inliers = 1.0;
    double n_pts = 1.0;
    double h = 0.0, k = 0.0, prob_reject_good_model = 0.0, log_eta = 0.0;
    int total_pts = u->data_size_;
    double new_eps = inl_nr * 1.0/total_pts;
    int stop_nr;

    for (int i = 0; i < u->min_sample_size; ++i) {
      n_inliers *= inl_nr - i;
      n_pts *= total_pts - i;
    }
    double prob_good_model = n_inliers/n_pts;

    if (prob_good_model < EPS) {
      stop_nr = u->max_hypotheses;
      goto out;
    } else if (1 - prob_good_model < EPS) {
      stop_nr = 1;
      goto out;
    }

    for (auto &i: wald_test_history_) {
      k += i.k_;
      h = compute_exp_sprt(new_eps, i.epsilon_, i.delta_);
      prob_reject_good_model = 1/(exp( h*log(i.A_) ));
      log_eta += (double) i.k_ * log( 1 - prob_good_model*(1-prob_reject_good_model) );
    }

    stop_nr = (int) ceil(k + ( log(1 - u->conf_th) - log_eta ) / log( 1-prob_good_model * (1-(1/decision_threshold_))));

  out:
    *stopping_nr = stop_nr;
    need_update_stopping_ = false;
  }

  int Losac::locally_optimize_soln(USAC *u, int best_inliers)
  {
    if (best_inliers < inner_sample_size / USAC_LOSAC_INNER_SAMPLE_SIZE_RATIO)
      return 0;

    int lo_sample_size = std::min(inner_sample_size, (int)(best_inliers * USAC_LOSAC_INNER_SAMPLE_SIZE_RATIO));
    vector<int> sample(lo_sample_size);
    vector<int> orig_inliers;
    vector<int> iter_inliers;

    // find all inliers less than threshold
    int lo_inliers = best_inliers;
    find_inliers(u->best_errs_, u->inlier_th_, &orig_inliers);

    u->results_.num_local_optimizations_++;

    vector<double> weights;
    double threshold_step_size = (threshold_multiplier * u->inlier_th_ - u->inlier_th_) / num_iterative_steps;

    for (int i = 0; i < inner_ransac_reps; ++i) {
      // generate non-minimal sample model and find inliers
      gen_uniform_rand_sample(orig_inliers.size(), lo_sample_size, &sample);
      
      for (int j = 0; j < lo_sample_size; ++j)
	sample[j] = orig_inliers[sample[j]];
      
      if (!u->gen_refined_model(sample))
	continue;
      if (!u->evaluate_model(0, &u->tmp_errs_, NULL, NULL))
	continue;
      find_inliers(u->tmp_errs_, threshold_multiplier*u->inlier_th_, &iter_inliers);
      if (!u->gen_refined_model(iter_inliers))
	continue;

      // iterative (reweighted) refinement - reduce threshold in steps, find new inliers and refit fundamental matrix
      // using weighted least-squares
      for (int j = 0; j < num_iterative_steps; ++j) {
	if (!u->evaluate_model(0, &u->tmp_errs_, NULL, NULL))
	  continue;
	find_inliers(u->tmp_errs_, (threshold_multiplier*u->inlier_th_) - (j+1)*threshold_step_size, &iter_inliers);
	u->find_weights(0, iter_inliers, &weights);
      
	if (!u->gen_refined_model(iter_inliers, &weights))
	  continue;
      }

      int inl_nr;
      if (!u->evaluate_model(0, &u->tmp_errs_, &inl_nr, NULL))
	continue;

      if (inl_nr > lo_inliers) {
	std::cout << "----Updated from " << lo_inliers << " to " << inl_nr <<"\n";
	lo_inliers = inl_nr;
	u->store_soln(0, lo_inliers);
      }
    }
    return lo_inliers;
  }

  void Losac::optimize_and_update(USAC *u)
  {
    int lo_inlier_count = locally_optimize_soln(u, u->results_.best_inlier_nr_);
    if (lo_inlier_count > u->results_.best_inlier_nr_ && lo_inlier_count > u->min_inls) {
      std::cout << "(" << u->results_.hyp_nr_ << ") Performing LO. Inlier count before: " << u->results_.best_inlier_nr_ << ", inlier count after: " << lo_inlier_count << std::endl;
      u->results_.best_inlier_nr_ = lo_inlier_count;
    }

    if (num_prev_best_inliers_ < u->results_.best_inlier_nr_) {
      // backup the old set of inliers for reference in LO
      num_prev_best_inliers_ = u->results_.best_inlier_nr_;
      prev_best_inliers_ = u->results_.inlier_flags_;
    }
  }

  std::vector<unsigned int> USAC::allSamples;
  
  void USAC::writeSamplesToFile(const std::string& fileName){
    std::ofstream ofs(fileName.c_str());
    for (int i=0; i<allSamples.size(); i++)
      ofs << allSamples[i]<<"\t";  
    ofs.close();
  }
  
  void USAC::readSamplesFromFile(const std::string& fileName){
    std::ifstream ifs(fileName.c_str());
    allSamples.clear();
    while (!ifs.eof()){
      unsigned int spl;
      ifs >> spl;
      if (!ifs.eof())
	allSamples.push_back(spl);
    }
    ifs.close();      
  }
  
  void USAC::init_params(const USACParams& p)
  {
    min_sample_size = p.min_sample_size;
    CHECK(min_sample_size > 0);
    resize_fill_vec(min_sample_, min_sample_size, -1);
    min_inls = p.min_inls;

    conf_th = p.conf_th;
    max_hypotheses = p.max_hypotheses;
    min_hypotheses = p.min_hypotheses;
    max_solns_per_sample = p.max_solns_per_sample;
    need_prevalidate_sample = p.need_prevalidate_sample;
    need_prevalidate_model = p.need_prevalidate_model;
    need_test_degeneracy = p.need_test_degeneracy;
    sampling_method = p.sampling_method;
    verif_method = p.verif_method;
    need_local_optim = p.need_local_optim;
    mle_sac = p.mle_sac;

  
    if (sampling_method == SAMP_PROSAC)
      prosac = std::make_unique<Prosac>(p.prosac_max_samples, p.prosac_beta, p.prosac_non_rand_conf, p.prosac_min_stop_len);

    if (verif_method == VERIF_SPRT) {
      sprt = std::make_unique<Sprt>(p.sprt_tM, p.sprt_mS, p.sprt_delta, p.sprt_epsilon);
      evaluate_model = std::bind(&Sprt::evaluate_model, sprt.get(), this, _1, _2, _3, _4);
    } else if (verif_method == VERIF_STANDARD) {
      evaluate_model = std::bind(&USAC::evaluate_model_ransac, this, _1, _2, _3, _4);
    }

    if (need_local_optim)
      losac = std::make_unique<Losac>(p.losac_inner_sample_size, p.losac_inner_ransac_reps, p.losac_threshold_multiplier, p.losac_num_iterative_steps);
  }

  void USAC::init_data(int data_size, double inlier_th, const vector<int> *prosac_sorted_indices)
  {
    results_.reset();
    allSamples.clear();
    currentSampleIdx = 0; 
  
    data_size_ = data_size;
    inlier_th_ = inlier_th;

    if (sampling_method == SAMP_PROSAC && prosac_sorted_indices == NULL) {
      sampling_method = SAMP_UNIFORM;
      prosac = NULL;
    }

    if (sampling_method == SAMP_PROSAC)
      prosac->init_data(data_size, min_sample_size, max_hypotheses, *prosac_sorted_indices);

    if (verif_method == VERIF_SPRT)
      sprt->init_data();

    if (need_local_optim)
      losac->init_data(data_size_);

    resize_fill_vec(tmp_errs_, data_size_);
    resize_fill_vec(best_errs_, data_size_);

    results_.resize(data_size_, min_sample_size, need_test_degeneracy);
    adaptive_stopping_nr_ = max_hypotheses;
  }

  IterState USAC::iterate(int iter_nr) {
    int min_sample_size = min_sample_.size();
    if (data_size_ < min_sample_size
	|| (sampling_method == SAMP_PROSAC &&
	    data_size_ < prosac->min_stop_len))
      return ITER_NOT_ENOUGH_DATA;

    IterState ret = ITER_CONTINUE;

    for (int cur_iter = 0; cur_iter < iter_nr; cur_iter++) {
    
      if ( (results_.hyp_nr_ >= fmin(adaptive_stopping_nr_, max_hypotheses)
	    && results_.hyp_nr_> min_hypotheses) )	
	return ITER_END;
      
      //std::cout << adaptive_stopping_nr_ << "\t" << results_.hyp_nr_ <<"\n";
      
      ++results_.hyp_nr_;
    
      if (sampling_method == SAMP_UNIFORM){
      
	gen_uniform_rand_sample(data_size_, min_sample_size, &min_sample_);
	//Save sample to allSample vector
	for (int i=0; i<min_sample_size; i++)
	  allSamples.push_back(min_sample_[i]);
      }
      else if (sampling_method == SAMP_PROSAC){
	prosac->gen_min_sample(this, &min_sample_);
      }
      else if (sampling_method == SAMP_FILE){
	if (currentSampleIdx >= allSamples.size())
	  return ITER_END;
	
	get_sample_from_list(allSamples, currentSampleIdx, min_sample_size, &min_sample_);
	
      }

      if (need_prevalidate_sample) {
	if (!validate_sample(min_sample_)) {
	  ++results_.rejected_sample_nr_;
	  continue;
	}
      }

      int num_solns = gen_min_sample_models(min_sample_);
      results_.model_nr_ += num_solns;

      bool update_best = false;
      for (int i = 0; i < num_solns; ++i) {
	if (need_prevalidate_model) {
	  if (!validate_model(i, min_sample_)) {
	    ++results_.rejected_model_nr_;
	    continue;
	  }
	}
	
	int inlier_count, test_nr;
	bool good = evaluate_model(i, &tmp_errs_, &inlier_count, &test_nr);
	results_.total_points_verified_ += test_nr;

	//If MLE_SAC is required, compute maximum likelyhood
	
	if (mle_sac){
	  double dSigma = inlier_th_/2;
	  double penalty = penaltyMLE(&tmp_errs_, dSigma);
	  //std::cout << penalty <<"\t" << inlier_count <<"\n";
	  if (verify_and_store_mlesac(i, penalty, inlier_count))
	    update_best = true;
	}
       
	
	if( !mle_sac  && ((verif_method == VERIF_STANDARD && verify_and_store(i, inlier_count, test_nr, good)) || (verif_method == VERIF_SPRT && sprt->verify_and_store(this, i, inlier_count, test_nr, good))))
	  update_best = true;
      }

      bool degenerate_model = false, upgrade_model = false;

      if (need_test_degeneracy && update_best) {
	test_soln_degeneracy(&degenerate_model, &upgrade_model);
	if (degenerate_model && upgrade_model) {
	  int upgrade_inliers = upgrade_degenerate_model(&tmp_errs_);
	  if (upgrade_inliers > results_.best_inlier_nr_ && upgrade_inliers > min_inls)
	    results_.best_inlier_nr_ = upgrade_inliers;
	}
      }

      if (need_local_optim && update_best){
	std::cout << "Running Local Update (LORANSAC)....\n";
	losac->optimize_and_update(this);
      }

      if (update_best) {
	if (sampling_method == SAMP_PROSAC){
	  adaptive_stopping_nr_ = prosac->update_stopping(this, results_.hyp_nr_);
	}
	else{
	  adaptive_stopping_nr_ = update_stopping(results_.best_inlier_nr_, data_size_, min_sample_size);
	  std::cout << " ----------Update stopping ------------- \n";
	  std::cout << " Best inls : " << results_.best_inlier_nr_ << " Data size : " << data_size_ << " Min sample size : " << min_sample_size << "\n";
	  std::cout  << "Iter nr : " << iter_nr << " N Stop : " << adaptive_stopping_nr_ << "\n";	  
	}
      }

      if (verif_method == VERIF_SPRT && sampling_method != SAMP_PROSAC)
	sprt->update_stopping(this, results_.best_inlier_nr_, &adaptive_stopping_nr_);
    }
    
    return ret;
  }

  bool USAC::solve()
  {
    clock_t tick = clock();
    auto flag = iterate(max_hypotheses);
    if (flag == ITER_NOT_ENOUGH_DATA)
      return false;
    CHECK(flag == ITER_END);
    double total_runtime_ = (clock() - tick) * 1.0 / CLOCKS_PER_SEC;
    std::cout << "Time: " << total_runtime_ << std::endl << std::endl;
    return true;
  }

  int USAC::update_stopping(int inl_nr, int total_pts, int sample_size)
  {
    double n_inliers = 1.0;
    double n_pts = 1.0;

    for (int i = 0; i < sample_size; ++i) {
      n_inliers *= inl_nr - i;
      n_pts *= total_pts - i;
    }
    double prob_good_model = n_inliers/n_pts;

    if ( prob_good_model < std::numeric_limits<double>::epsilon() )
      return max_hypotheses;
    else if ( 1 - prob_good_model < std::numeric_limits<double>::epsilon() )
      return 1;
    else{
      int maxHyp = (int) ceil(log(1-conf_th)/log(1-prob_good_model));
      if (maxHyp < 0) return max_hypotheses;      
      return maxHyp;
    }
  }


  void USAC::store_soln(int model_idx, int inl_nr)
  {
    results_.best_inlier_nr_ = inl_nr;
    for (int i = 0; i < (int) tmp_errs_.size(); ++i) {
      if (tmp_errs_[i] < inlier_th_)
	results_.inlier_flags_[i] = true;
      else
	results_.inlier_flags_[i] = false;
    }

    results_.best_sample_ = min_sample_;
    best_errs_ = tmp_errs_;
    store_model(model_idx, inl_nr);
  }

  bool USAC::evaluate_model_ransac(int model_idx, vector<double> *errs, int *inlier_nr, int *test_nr)
  {
    int inl_nr = 0;
    for (int i = 0; i < (int) errs->size(); ++i) {
      double err = compute_error(model_idx, i);
      (*errs)[i] = err;
      if (err < inlier_th_)
	inl_nr++;
    }
    if (inlier_nr)
      *inlier_nr = inl_nr;
    if (test_nr)
      *test_nr = errs->size();
    return true;
  }

  bool USAC::verify_and_store(int model_idx, int inlier_nr, int test_nr, bool good)
  {
    if (inlier_nr > results_.best_inlier_nr_ && inlier_nr > min_inls) {
      results_.best_inlier_nr_ = inlier_nr;
      store_soln(model_idx, results_.best_inlier_nr_);
      return true;
    }
    return false;
  }

  bool USAC::verify_and_store_mlesac(int model_idx, double penalty, int inlier_nr){
    if (penalty < results_.min_mle_penalty_){
      results_.min_mle_penalty_ = penalty;
      results_.best_inlier_nr_ = inlier_nr;
      store_soln(model_idx, results_.best_inlier_nr_);
      return true;      
    }
    return false;
  }
								    

  double USAC::penaltyMLE(vector<double> *errs, double dSigma){
    
    double dMix = 0.5;
    double dResidInlierProb, dResidOutlierProb;       
    //    double GlobMaxResid = 500.0;
    double GlobMaxResid = 1.0;

    double dInlierProb =  0.0;
    for (int EM_iter=0; EM_iter<nEMIter; EM_iter++){
      double dInlierProb = 0.0;
      for (int i=0; i<errs->size(); i++){
	dResidInlierProb = dMix * exp(-(*errs)[i] * (*errs)[i]/ (2*dSigma*dSigma))/(dSigma*sqrt(2*M_PI));
	dResidOutlierProb = (1.0-dMix)/GlobMaxResid;
	dInlierProb += dResidInlierProb/(dResidInlierProb + dResidOutlierProb);
      }
      dMix = dInlierProb/errs->size();      
    }
    
    /*Find Log-likelyhood*/
    double dCurPenalty = 0;
    for (int i=0; i<errs->size(); i++){
      dResidInlierProb = dMix * exp(-(*errs)[i] * (*errs)[i]/ (2*dSigma*dSigma))/(dSigma*sqrt(2*M_PI));
      dResidOutlierProb = (1.0-dMix)/GlobMaxResid;
      
      dCurPenalty += -log(dResidInlierProb + dResidOutlierProb);	
    }
    return dCurPenalty;
  }

  
  
} // namespace usac
