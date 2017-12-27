#include "AffineEstimator.hh"
#include <iostream>



using namespace std;

namespace usac{
  
  bool AffineEstimator::init_problem(Eigen::MatrixXd *A, Eigen::MatrixXd *b, Eigen::MatrixXd *c, Eigen::VectorXd *d){
    AA = A; bb = b; cc = c; dd = d;
    npts = bb->cols();
    final_model_params_.resize(6);
    resize_fill_vec(models_, max_solns_per_sample);
    
    models_[0] = new double[6];
    return true;        
  }
  
  int AffineEstimator::gen_min_sample_models(const vector<int>& sample){

    int sampleSize = sample.size();
    Eigen::MatrixXd xx(sampleSize*2, 6);
    Eigen::VectorXd yy(sampleSize*2);

    for (int i=0; i<sampleSize; i++){
      int spl = sample[i];
      xx.row(2*i) = -(*AA).row(2*spl);
      xx.row(2*i+1)= -(*AA).row(2*spl);
      yy(2*i) = (*bb)(0,spl);
      yy(2*i+1) = (*bb)(1,spl);      
    }
    
    Eigen::VectorXd theta = xx.colPivHouseholderQr().solve(yy);
    /*Pass solution to model*/
    for (int i=0; i<6; i++){
      *(models_[0]+i)=theta[i];
    }
    return 1;      
    
  }

  double AffineEstimator::compute_error(int model_idx, int sample_idx){
    double *model = models_[model_idx];
    double sum = 0;
    Eigen::VectorXd tt(6);
    for (int i=0; i<6; i++){
      tt(i) = *(models_[0]+i);
    }
    sum = fabs( (*AA).row(2*sample_idx)*tt + (*bb)(0,sample_idx)) +
      fabs((*AA).row(2*sample_idx+1)*tt + (*bb)(1, sample_idx));    
    return sum;
  }

  void AffineEstimator::store_model(int model_idx, int inlier_nr){
    for (int i=0; i<6; i++)
      final_model_params_[i] = *(models_[model_idx]+i);
  }

  bool AffineEstimator::gen_refined_model(const vector<int>&sample, const vector<double>* weights){
    int larger_sample_size = sample.size();
    
    if (larger_sample_size == 0)
      return false;
    
    Eigen::MatrixXd xx(2*larger_sample_size, 6);
    Eigen::VectorXd yy(2*larger_sample_size);
    Eigen::MatrixXd W(2*larger_sample_size, 2*larger_sample_size);

    for (int i=0; i<sample.size(); i++){
      int spl = sample[i];
      xx.row(2*i) = (*AA).row(2*spl);
      xx.row(2*i+1) = (*AA).row(2*spl+1);
      yy[2*i] = (*bb)(0, spl);
      yy[2*i+1] = (*bb)(1,spl);
      double wi = 1;
      if (weights !=NULL){
	wi = (*weights)[i];
      }
      W(2*i, 2*i) = wi;
      W(2*i+1, 2*i+1)=wi;      
    }
      
    Eigen::MatrixXd Xp = W*xx;
    Eigen::VectorXd Yp = W*yy;
    Eigen::VectorXd theta = Xp.jacobiSvd(ComputeThinU | ComputeThinV).solve(Yp);

    /*Pass solution to Model*/
    for (int i=0; i<6; i++){
      *(models_[0]+i) = theta[i];
    }    
  }

  void AffineEstimator::find_weights(int model_idx, const vector<int>&inliers, vector<double>*weights){
    int inlier_nr = inliers.size();
    (*weights).resize(inlier_nr);
    Eigen::VectorXd tt(6);
    for (int i=0; i<6; i++){
      tt[i] = *(models_[model_idx]+i);
    }
        
    for (int i=0; i < inlier_nr; i++){
      int spl = inliers[i];
      double r = fabs( (*AA).row(2*spl)*tt + (*bb)(0,spl)) + fabs( (*AA).row(2*spl+1)*tt + (*bb)(0, spl));
      (*weights)[i] = 1.0/r;
    }  
  }

  bool AffineEstimator::validate_model(int modelIndex, const vector<int>& sample){
    return true;
  }
  
  
  
  
  
}








