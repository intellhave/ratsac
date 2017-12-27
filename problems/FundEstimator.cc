#include <iostream>
#include <fstream>
#include <string>
#include "MathFunctions.hh"
#include "FundFunctions.hh"
#include "HomographyFunctions.hh"
#include "FundEstimator.hh"


namespace usac {

  bool FundEstimator::init_problem(const FundParams &para, int nr, double* pointData)
  {
    input_points_denorm_ = pointData;
    input_points_       = new double[6 * nr];
    if (input_points_denorm_ == NULL) {
      std::cerr << "Input point data not properly initialized" << std::endl;
      return false;
    }
    if (input_points_ == NULL)
      {
	std::cerr << "Could not allocate storage for normalized data points" << std::endl;
	return false;
      }
    
    FTools::normalizePoints(input_points_denorm_, input_points_, nr, m_T1_, m_T2_);
    MathTools::mattr(m_T2_trans_, m_T2_, 3, 3);

    final_model_params_.clear(); final_model_params_.resize(9);
    resize_fill_vec(models_, max_solns_per_sample);
    resize_fill_vec(models_denorm_, max_solns_per_sample);
    for (int i = 0; i < max_solns_per_sample; ++i) {
      models_[i] = new double[9];
      models_denorm_[i] = new double[9];
    }

    data_matrix_ = new double[9*nr];
    FTools::computeDataMatrix(data_matrix_, nr, input_points_);

    resize_fill_vec(degen_outlier_flags_, nr, 0);
    if (need_test_degeneracy) {
      degen_homog_threshold_ = para.h_degen_th * para.h_degen_th;
      degen_max_upgrade_samples_ = para.max_upgrade_samples;
      resize_fill_vec(degen_final_model_params_, 9);
      degen_data_matrix_ = new double[2*9*nr];
      HTools::computeDataMatrix(degen_data_matrix_, nr, input_points_);
    }
    matrix_decomposition_method = para.decomp_alg;
    return true;
  }

  void FundEstimator::cleanup_problem()
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
    degen_outlier_flags_.clear();
  }

  int FundEstimator::gen_min_sample_models(const vector<int> &sample)
  {
    double A[9*9];
    int nsols = 0;

    // form the matrix of equations for this minimal sample
    double *src_ptr;
    double *dst_ptr = A;
    for (auto i: sample) {
      src_ptr = data_matrix_ + i;
      for (int j = 0; j < 9; ++j)
	{
	  *dst_ptr = *src_ptr;
	  ++dst_ptr;
	  src_ptr += data_size_;
	}
    }

    // LU/QR factorization
    double sol[9*9];
    double poly[4], roots[3];
    double *f1, *f2;
    int nullbuff [18];
    f1 = sol;
    f2 = sol+9;
    if (matrix_decomposition_method == DECOMP_QR) {
      FTools::nullspaceQR7x9(A, sol);

    } else if (matrix_decomposition_method == DECOMP_LU) {
      for (int i = 7*9; i < 9*9; ++i)
	A[i] = 0.0;

      int nullsize = FTools::nullspace(A, f1, 9, nullbuff);
      if (nullsize != 2)
	return 0;
    }

    // solve polynomial
    FTools::makePolynomial(f1, f2, poly);
    nsols = FTools::rroots3(poly, roots);

    // form up to three fundamental matrices
    double T2_F[9];
    for (int i = 0; i < nsols; ++i) {
      for (int j = 0; j < 9; ++j)
	*(models_[i]+j) = f1[j] * roots[i] + f2[j] * (1 -roots[i]);

      // store denormalized version as well
      MathTools::mmul(T2_F, m_T2_trans_, models_[i], 3);
      MathTools::mmul(models_denorm_[i], T2_F, m_T1_, 3);
    }

    return nsols;
  }


  bool FundEstimator::gen_refined_model(const vector<int>& sample, const vector<double> * weights)
  {
    // form the matrix of equations for this non-minimal sample
    int numPoints = sample.size();
    double *A = new double[numPoints*9];
    double *src_ptr;
    double *dst_ptr = A;
    for (int i = 0; i < numPoints; ++i)
      {
	src_ptr = data_matrix_ + sample[i];
	for (int j = 0; j < 9; ++j)
	  {
	    if (weights == NULL)
	      *dst_ptr = *src_ptr;
	    else
	      *dst_ptr = (*src_ptr)*(*weights)[i];

	    ++dst_ptr;
	    src_ptr += data_size_;
	  }
      }

    double Cv[9*9];
    FTools::formCovMat(Cv, A, numPoints, 9);

    double V[9*9], D[9], *p;
    MathTools::svdu1v(D, Cv, 9, V, 9);

    int j = 0;
    for (int i = 1; i < 9; ++i)
      {
	if (D[i] < D[j])
	  {
	    j = i;
	  }
      }
    p = V + j;

    for (int i = 0; i < 9; ++i)
      {
	*(models_[0]+i) = *p;
	p += 9;
      }
    FTools::singulF(models_[0]);
    // store denormalized version as well
    double T2_F[9];
    MathTools::mmul(T2_F, m_T2_trans_, models_[0], 3);
    MathTools::mmul(models_denorm_[0], T2_F, m_T1_, 3);

    delete[] A;

    return true;
  }

  bool FundEstimator::validate_model(int modelIndex, const vector<int> &sample)
  {
    // check oriented constraints
    double e[3], sig1, sig2;
    FTools::computeEpipole(e, models_[modelIndex]);

    sig1 = FTools::getOriSign(models_[modelIndex], e, input_points_ + 6*sample[0]);
    for(int i = 1; i < (int)min_sample_.size(); ++i) {
      sig2 = FTools::getOriSign(models_[modelIndex], e, input_points_ + 6*sample[i]);
      if (sig1 * sig2 < 0)
	return false;
    }
    return true;
  }

  double FundEstimator::compute_error(int model_idx, int sample_idx)
  {
    double *model = models_denorm_[model_idx];
    double *pt = input_points_denorm_ + 6*sample_idx;
    double rxc = (*model) * (*(pt+3)) + (*(model+3)) * (*(pt+4)) + (*(model+6));
    double ryc = (*(model+1)) * (*(pt+3)) + (*(model+4)) * (*(pt+4)) + (*(model+7));
    double rwc = (*(model+2)) * (*(pt+3)) + (*(model+5)) * (*(pt+4)) + (*(model+8));
    double r =((*(pt)) * rxc + (*(pt+1)) * ryc + rwc);
    double rx = (*model) * (*(pt)) + (*(model+1)) * (*(pt+1)) + (*(model+2));
    double ry = (*(model+3)) * (*(pt)) + (*(model+4)) * (*(pt+1)) + (*(model+5));
    return r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry);

  }

  void FundEstimator::test_soln_degeneracy(bool* degenerateModel, bool* upgradeModel)
  {
    *degenerateModel = false;
    *upgradeModel = false;

    // make up the tuples to be used to check for degeneracy
    int degen_sample_indices[] = {0, 1, 2, 3,
				  3, 4, 5, 6,
				  0, 1, 5, 6,
				  0, 2, 4, 5,
				  1, 2, 4, 6,
				  0, 3, 4, 6,
				  1, 3, 4, 5,
				  2, 3, 5, 6};

    // the above tuples need to be tested on the remaining points for each case
    int test_point_indices[] = {4, 5, 6,
				0, 1, 2,
				2, 3, 4,
				1, 3, 6,
				0, 3, 5,
				1, 2, 5,
				0, 2, 6,
				0, 1, 4};

    int *sample_pos = degen_sample_indices;
    int *test_pos = test_point_indices;
    double h[9];
    double T2_inv[9], T2_H[9];
    for (int i = 0; i < 9; ++i)
      T2_inv[i] = m_T2_[i];

    MathTools::minv(T2_inv, 3);

    vector<int> sample(7), test(3);
    vector<double> errs;
    for(int i = 0; i < 8; ++i)
      {
	// compute H from the current set of 4 points
	for (int j = 0; j < 4; ++j)
	  sample[j] = min_sample_[sample_pos[j]];

	FTools::computeHFromMinCorrs(sample, 4, data_size_, degen_data_matrix_, h);
	MathTools::mmul(T2_H, T2_inv, h, 3);
	MathTools::mmul(h, T2_H, m_T1_, 3);

	// check test points to see how many are consistent
	for (int j = 0; j < 3; ++j)
	  test[j] = min_sample_[test_pos[j]];

	int num_inliers = FTools::getHError(errs, input_points_denorm_, data_size_, h, degen_homog_threshold_, &test);
	for (int j = 0, count = 4; j < 3; ++j)
	  if (errs[j] < degen_homog_threshold_)
	    sample[count++] = test[j];

	// if at least 1 inlier in the test points, then h-degenerate sample found
	if (num_inliers > 0)
	  {
	    // find inliers from all data points
	    num_inliers = FTools::getHError(errs, input_points_denorm_, data_size_, h, degen_homog_threshold_);
	    //std::cout << "Degenerate sample found with " << num_inliers << " inliers" << std::endl;

	    // refine with least squares fit
	    int count = 0;
	    vector<int> inlier_sample(num_inliers);
	    for (int j = 0; j < data_size_; ++j)
	      if (errs[j] < degen_homog_threshold_)
		inlier_sample[count++] = j;

	    FTools::computeHFromCorrs(inlier_sample, inlier_sample.size(), data_size_, degen_data_matrix_, h);
	    MathTools::mmul(T2_H, T2_inv, h, 3);
	    MathTools::mmul(h, T2_H, m_T1_, 3);

	    // find support of homography
	    num_inliers = FTools::getHError(errs, input_points_denorm_, data_size_, h, degen_homog_threshold_);
	    //std::cout << "Refined model has " << num_inliers << " inliers" << std::endl;
#if 1
	    if (num_inliers < results_.best_inlier_nr_/5)
	      {
		sample_pos += 4;
		test_pos += 3;
		continue;
	      }
#endif
	    // set flag
	    *degenerateModel = true;

	    // if largest degenerate model found so far, store results
	    if (num_inliers > results_.degen_inlier_nr_)
	      {
		// set flag
		std::cerr << "upgrade model " << std::endl;
		*upgradeModel = true;

		// refine with final least squares fit
		count = 0;
		inlier_sample.resize(num_inliers);
		for (int j = 0; j < data_size_; ++j)
		  if (errs[j] < degen_homog_threshold_)
		    inlier_sample[count++] = j;

		FTools::computeHFromCorrs(inlier_sample, inlier_sample.size(), data_size_, degen_data_matrix_, h);

		results_.degen_inlier_nr_ = num_inliers;
		// store homography
		for (int j = 0; j < 9; ++j)
		  degen_final_model_params_[j] = h[j];

		// store inliers and outliers - for use in model completion
		for (int j = 0; j < data_size_; ++j)
		  {
		    if (errs[j] < degen_homog_threshold_)
		      {
			results_.degen_inlier_flags_[j] = 1;
			degen_outlier_flags_[j] = 0;
		      }
		    else
		      {
			degen_outlier_flags_[j] = 1;
			results_.degen_inlier_flags_[j] = 0;
		      }
		  }
		// store the degenerate points from the minimal sample
		results_.degen_sample_ = sample;

	      } // end store denerate results
	  } // end check for one model degeneracy

	sample_pos += 4;
	test_pos += 3;

      } // end check for all combinations in the minimal sample
  }


  int FundEstimator::upgrade_degenerate_model(vector<double> *errs)
  {
    int best_upgrade_inliers = results_.best_inlier_nr_;
    int num_outliers = data_size_ - results_.degen_inlier_nr_;

    if (num_outliers < 2)
      return 0;

    vector<int> outlier_indices(num_outliers);
    int count = 0;
    for (int i = 0; i < data_size_; ++i)
      if (degen_outlier_flags_[i])
	outlier_indices[count++] = i;

    vector<int> outlier_sample(2);

    double* pt1_index, *pt2_index;
    double x1[3], x1p[3], x2[3], x2p[3];
    double temp[3], l1[3], l2[3], ep[3];
    double skew_sym_ep[9];
    double T2_F[9];
    std::cerr  << "num outliers "<< num_outliers << std::endl;
    for (int i = 0; i < degen_max_upgrade_samples_; ++i) {
      gen_uniform_rand_sample(num_outliers, 2, &outlier_sample);

      pt1_index = input_points_ + 6*outlier_indices[outlier_sample[0]];
      pt2_index = input_points_ + 6*outlier_indices[outlier_sample[1]];

      x1[0]  = pt1_index[0]; x1[1]  = pt1_index[1]; x1[2]  = 1.0;
      x1p[0] = pt1_index[3]; x1p[1] = pt1_index[4]; x1p[2] = 1.0;
      x2[0]  = pt2_index[0]; x2[1]  = pt2_index[1]; x2[2]  = 1.0;
      x2p[0] = pt2_index[3]; x2p[1] = pt2_index[4]; x2p[2] = 1.0;

      MathTools::vmul(temp, &degen_final_model_params_[0], x1, 3);
      MathTools::crossprod(l1, temp, x1p, 1);

      MathTools::vmul(temp, &degen_final_model_params_[0], x2, 3);
      MathTools::crossprod(l2, temp, x2p, 1);

      MathTools::crossprod(ep, l1, l2, 1);

      MathTools::skew_sym(skew_sym_ep, ep);
      MathTools::mmul(models_[0], skew_sym_ep, &degen_final_model_params_[0], 3);
      MathTools::mmul(T2_F, m_T2_trans_, models_[0], 3);
      MathTools::mmul(models_denorm_[0], T2_F, m_T1_, 3);

      int num_inliers;
      evaluate_model(0, errs, &num_inliers, NULL);

      if (num_inliers > best_upgrade_inliers)
	{
	  std::cerr << "upgrade inliers " << num_inliers << std::endl;
	  results_.degen_sample_[5] = outlier_indices[outlier_sample[0]];
	  results_.degen_sample_[6] = outlier_indices[outlier_sample[1]];
	  min_sample_ = results_.degen_sample_;
	  store_soln(0, num_inliers);
	  best_upgrade_inliers = num_inliers;

	  int count = 0;
	  for (size_t j = 0; j < outlier_indices.size(); ++j)
	    if ((*errs)[outlier_indices[j]] < inlier_th_)
	      ++count;

	  int num_samples = update_stopping(count, num_outliers, 2);
	  if (num_samples < degen_max_upgrade_samples_)
	    degen_max_upgrade_samples_ = num_samples;
	}
    }

    std::cout << "Upgraded model has " << best_upgrade_inliers << " inliers" << std::endl;
    return best_upgrade_inliers;
  }



  //Find weight overwirte
  void FundEstimator::find_weights(int model_idx, const vector<int>&inliers, vector<double>*weights){
    double rx, ry, ryc, rxc;
    double* model = models_[model_idx];
    double* pt;
    int pt_index;

    int inlier_nr = inliers.size();
    for (int i = 0; i < inlier_nr; ++i)
      {
	// get index of point to be verified
	pt_index = inliers[i];

	// compute weight (ref: torr dissertation, eqn. 2.25)
	pt = input_points_ + 6*pt_index;
	rxc = (*model) * (*(pt+3)) + (*(model+3)) * (*(pt+4)) + (*(model+6));
	ryc = (*(model+1)) * (*(pt+3)) + (*(model+4)) * (*(pt+4)) + (*(model+7));
	rx = (*model) * (*(pt)) + (*(model+1)) * (*(pt+1)) + (*(model+2));
	ry = (*(model+3)) * (*(pt)) + (*(model+4)) * (*(pt+1)) + (*(model+5));
	double sqrtd = sqrt(rxc*rxc + ryc*ryc + rx*rx + ry*ry);
	(*weights).push_back(1/sqrtd);
      }
    
  }
  
  //

  // ============================================================================================
  // find_weights: given model and points, compute weights to be used in local optimization
  // ============================================================================================
  void FundEstimator::find_weights(int model_idx, const vector<int>& inliers,
				   int inlier_nr, double* weights)
  {
    double rx, ry, ryc, rxc;
    double* model = models_[model_idx];
    double* pt;
    int pt_index;

    for (int i = 0; i < inlier_nr; ++i)
      {
	// get index of point to be verified
	pt_index = inliers[i];

	// compute weight (ref: torr dissertation, eqn. 2.25)
	pt = input_points_ + 6*pt_index;
	rxc = (*model) * (*(pt+3)) + (*(model+3)) * (*(pt+4)) + (*(model+6));
	ryc = (*(model+1)) * (*(pt+3)) + (*(model+4)) * (*(pt+4)) + (*(model+7));
	rx = (*model) * (*(pt)) + (*(model+1)) * (*(pt+1)) + (*(model+2));
	ry = (*(model+3)) * (*(pt)) + (*(model+4)) * (*(pt+1)) + (*(model+5));

	weights[i] = 1/sqrt(rxc*rxc + ryc*ryc + rx*rx + ry*ry);
      }
  }

  void FundEstimator::store_model(int model_idx, int inlier_nr)
  {
    for (int i = 0; i < 9; ++i)
      final_model_params_[i] = *(models_denorm_[model_idx]+i);
  }

} //namespace usac
