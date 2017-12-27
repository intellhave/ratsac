#include <iostream>
#include <fstream>
#include <string>
#include "LinearEstimator.hh"
#include <armadillo>
#include <set>
#include "Chebyshev.hh"
using namespace arma;
using namespace std;

namespace usac{

  bool LinearEstimator::init_problem(const LinearParams &para, Eigen::MatrixXd* x_, Eigen::VectorXd* y_)
  {
    xData = x_;
    yData = y_;
    npts = (*xData).rows();
    d = (*xData).cols();
    
    //final_model_params_.clear();
    final_model_params_.resize(d);

    resize_fill_vec(models_, max_solns_per_sample);
    resize_fill_vec(models_denorm_, max_solns_per_sample);

    for (int i=0; i < max_solns_per_sample; i++){
      models_[i] = new double[d];
      models_denorm_[i] = new double[d];
    }
    return true;
    
  }
 
  int LinearEstimator::gen_min_sample_models(const vector<int> &sample){
    // Perform Chebyshev fitting here
    //((x'*x)^-1)*x'*(y)
    // Temporary do least square fit
    int sampleSize = sample.size();
    Eigen::MatrixXd xSampled(sampleSize, d);
    Eigen::VectorXd ySampled(sampleSize);
    
    for (int i=0; i<sample.size(); i++){      
      xSampled.row(i) = (*xData).row(sample[i]);
      ySampled[i] = (*yData)[sample[i]];
    }
    //vec theta = arma::solve(xSampled, ySampled);
    // vec theta0 = randu<vec>(d);
    // set<int> basics;
    // double minimax;
    
    // vec theta = Chebyshev::Descend(x, y, inlier_th, theta0, basics, minimax);

    //double inlier_th = 0.1;
    //Eigen::VectorXd theta0 = Eigen::VectorXd::Random(d);
    //double minimax;
    //set<int>solutionBasics;
    //Eigen::VectorXd theta = Chebyshev::Descend(xSampled, ySampled, inlier_th, theta0, solutionBasics, minimax);
    Eigen::VectorXd theta = xSampled.colPivHouseholderQr().solve(ySampled);
    
    
    //Pass solution to model
    for (int i=0; i<d; i++){
      *(models_[0]+i)=theta[i];
    }
    
    return 1;
    
    
  }

  double LinearEstimator::compute_error(int model_idx, int sample_idx){
    
    double *model = models_[model_idx];
    double sum = 0;
    for (int i=0; i<d; i++)
      sum = sum + (*xData)(sample_idx,i)*(*(model+i));
    
    double y1 = (*yData)[sample_idx];
    
    double r = sum - y1;
    if (r<0) r = -r;
    return r;    
  }

  void LinearEstimator::store_model(int model_idx, int inlier_nr){
    for (int i=0; i<d; ++i)
      final_model_params_[i] = *(models_[model_idx]+i);
    
  }


  bool LinearEstimator::gen_refined_model(const vector<int>&sample, const vector<double>*weights){
    //There is a larger-than-minimal sample together with the weights
    // Estimate the weighted least square fit;
       
    int larger_sample_size = sample.size();

    if (larger_sample_size==0)
      return false;
    
    Eigen::MatrixXd xSampled(larger_sample_size, d);
    Eigen::VectorXd ySampled(larger_sample_size);
    Eigen::MatrixXd W(larger_sample_size, larger_sample_size);
   
    //Fill the Data matrix and weight matrix
    for (int i=0; i<sample.size(); i++){
      //Data
      xSampled.row(i) = (*xData).row(sample[i]);
      ySampled[i] = (*yData)[sample[i]];
      //Weight Matrix
      double wi= 1;
      if (weights !=NULL){
	wi = (*weights)[i];
	//cout << wi << "  ";
      }
      W(i,i) = wi;
    }
    //Compute theta
    
    
    //Eigen::VectorXd theta = (xSampled.transpose()*W*xSampled).inverse()*xSampled.transpose()*W*ySampled;
    //Compute using SVD
    Eigen::MatrixXd Xp = W*xSampled;
    Eigen::VectorXd Yp = W*ySampled;
    Eigen::VectorXd theta = Xp.jacobiSvd(ComputeThinU | ComputeThinV).solve(Yp);
    //A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b)
    
    
    //Pass solution to model
    for (int i=0; i<d; i++){
      *(models_[0]+i) = theta[i];
    }
    
    return true;
  }

  bool LinearEstimator::validate_model(int modelIndex, const vector<int>&sample){
    return true;
  }


  void LinearEstimator::find_weights(int model_idx, const vector<int>&inliers,
				     vector<double>* weights){

    int inlier_nr = inliers.size();
    (*weights).resize(inlier_nr);
    
    Eigen::VectorXd theta(d);
    for (int i=0; i<d; i++) theta[i] = *(models_[model_idx]+i);
    for (int i=0; i<inlier_nr; i++){
      double r = ((*yData)[inliers[i]] -(*xData).row(inliers[i])*theta);
      (*weights)[i] = 1.0/abs(r);
    }
      

  }
  
  
  
}
