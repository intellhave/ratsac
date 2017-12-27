#ifndef HOMOGESTIMATOR_H
#define HOMOGESTIMATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include "MathFunctions.hh"
#include "FundFunctions.hh"
#include "HomographyFunctions.hh"
#include "HomogEstimator.hh"

namespace usac {

HomogEstimator::~HomogEstimator()
{
  cleanup_problem();
}

bool HomogEstimator::init_problem(int nr, double* pointData)
{
  // copy pointer to input data
  input_points_denorm_ = pointData;
  input_points_       = new double[6 * nr];
  if (input_points_denorm_ == NULL) {
    std::cerr << "Input point data not properly initialized" << std::endl;
    return false;
  }
  if (input_points_ == NULL) {
    std::cerr << "Could not allocate storage for normalized data points" << std::endl;
    return false;
  }

  // normalize input data
  // following this, input_points_ has the normalized points and input_points_denorm_ has
  // the original input points
  FTools::normalizePoints(input_points_denorm_, input_points_, nr, m_T1_, m_T2_);
  for (int i = 0; i < 9; ++i)
    m_T2inv_[i] = m_T2_[i];

  MathTools::minv(m_T2inv_, 3);

  // allocate storage for models
  final_model_params_.clear(); final_model_params_.resize(9);
  resize_fill_vec(models_, max_solns_per_sample);
  resize_fill_vec(models_denorm_, max_solns_per_sample);

  for (int i = 0; i < max_solns_per_sample; ++i) {
    models_[i] = new double[9];
    models_denorm_[i] = new double[9];
  }

  // precompute the data matrix
  data_matrix_ = new double[18*nr];  // 2 equations per correspondence
  HTools::computeDataMatrix(data_matrix_, nr, input_points_);

  return true;
}

void HomogEstimator::cleanup_problem()
{
  if (input_points_) { delete[] input_points_; input_points_ = NULL; }
  if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
  for (size_t i = 0; i < models_.size(); ++i)
    if (models_[i])
      delete[] models_[i];

  models_.clear();
  for (size_t i = 0; i < models_denorm_.size(); ++i)
    if (models_denorm_[i])
      delete[] models_denorm_[i];

  models_denorm_.clear();
}


int HomogEstimator::gen_min_sample_models(const vector<int> &sample)
{
  double A[8*9];
  double At[9*8];

  // form the matrix of equations for this minimal sample
  double *src_ptr;
  double *dst_ptr = A;
  for (auto i: sample) {
    for (int j = 0; j < 2; ++j) {
      src_ptr = data_matrix_ + 2*i + j;
      for (int k = 0; k < 9; ++k) {
        *dst_ptr = *src_ptr;
        ++dst_ptr;
        src_ptr += 2 * data_size_;
      }
    }
  }

  MathTools::mattr(At, A, 8, 9);

  double D[9], U[9*9], V[8*8], *p;
  MathTools::svduv(D, At, U, 9, V, 8);
  p = U + 8;

  double T2_H[9];
  for (int i = 0; i < 9; ++i)
  {
    *(models_[0]+i) = *p;
    p += 9;
  }
  MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
  MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3);

  for (int i = 0; i < 9; ++i)
    inv_model[i] = models_denorm_[0][i];
  MathTools::minv(inv_model, 3);
  return 1;
}


bool HomogEstimator::gen_refined_model(const vector<int>& sample, const vector<double>  *weights)
{
  CHECK(!weights || weights->size() == sample.size());
  int numPoints = sample.size();
  double *A = new double[numPoints*2*9];
  double *src_ptr;
  double *dst_ptr = A;
  for (int i = 0; i < numPoints; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      src_ptr = data_matrix_ + 2*sample[i] + j;
      for (int k = 0; k < 9; ++k)
      {
        if (weights == NULL)
          *dst_ptr = *src_ptr;
        else
          *dst_ptr = (*src_ptr) * (*weights)[i];

        ++dst_ptr;
        src_ptr += 2*data_size_;
      }
    }
  }

  // decompose
  double V[9*9], D[9], *p;
  MathTools::svdu1v(D, A, 2*numPoints, V, 9);

  int j = 0;
  for (int i = 1; i < 9; ++i)
  {
    if (D[i] < D[j])
      j = i;
  }
  p = V + j;

  for (int i = 0; i < 9; ++i)
  {
    *(models_[0]+i) = *p;
    p += 9;
  }
  double T2_H[9];
  MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
  MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3);

  for (int i = 0; i < 9; ++i)
    inv_model[i] = models_denorm_[0][i];
  MathTools::minv(inv_model, 3);

  delete[] A;

  return true;
}

bool HomogEstimator::validate_sample(const vector<int> &sample)
{
  // check oriented constraints
  double p[3], q[3];
  double *a, *b, *c, *d;

  a = input_points_ + 6*sample[0];
  b = input_points_ + 6*sample[1];
  c = input_points_ + 6*sample[2];
  d = input_points_ + 6*sample[3];

  HTools::crossprod(p, a, b, 1);
  HTools::crossprod(q, a+3, b+3, 1);

  if ((p[0]*c[0]+p[1]*c[1]+p[2]*c[2])*(q[0]*c[3]+q[1]*c[4]+q[2]*c[5])<0)
    return false;
  if ((p[0]*d[0]+p[1]*d[1]+p[2]*d[2])*(q[0]*d[3]+q[1]*d[4]+q[2]*d[5])<0)
    return false;

  HTools::crossprod(p, c, d, 1);
  HTools::crossprod(q, c+3, d+3, 1);

  if ((p[0]*a[0]+p[1]*a[1]+p[2]*a[2])*(q[0]*a[3]+q[1]*a[4]+q[2]*a[5])<0)
    return false;
  if ((p[0]*b[0]+p[1]*b[1]+p[2]*b[2])*(q[0]*b[3]+q[1]*b[4]+q[2]*b[5])<0)
    return false;

  return true;
}

double HomogEstimator::compute_error(int model_idx, int sample_idx)
{
  double h_x[3], h_inv_xp[3];
  double *model = models_denorm_[model_idx];
  double *pt = input_points_denorm_ + 6*sample_idx;
  MathTools::vmul(h_x, model, pt, 3);
  MathTools::vmul(h_inv_xp, inv_model, pt+3, 3);
  
  double err1 = 0.0, err2 = 0.0;
  double err = 0.0;
  for (int j = 0; j < 2; ++j) {
    err += fabs(h_x[j]/h_x[2] - pt[3+j]);
    err1 += (h_x[j]/h_x[2] - pt[3+j]) * (h_x[j]/h_x[2] - pt[3+j]);    
    err2 += (h_inv_xp[j]/h_inv_xp[2] - pt[j]) * (h_inv_xp[j]/h_inv_xp[2] - pt[j]);
  }

  Eigen::VectorXd theta(8);
  for (int i=0; i<8; i++)
    theta[i] = *(model+i)/(*(model+8));

  //std::cout << "CHECK ERROR FUNCITON HOMOG!!!!!!!!!!!!!!!!!!!!!\n";
  //return err1+err2;
  
  double err3 = QuasiconvexTools::computeError(*AA, *bb, *cc, *dd, theta, sample_idx);
  return err3;
  //std::cout<<"Err : " << err <<"\t" << err3 << "\n";
  
  //bool isinls = QuasiconvexTools::isInlier(*xData, *yData, theta, sample_idx);
  //cout << "Is inls " <<isinls << "\n";

  // if (err <= inlier_th_ && !isinls && err3>0){
  //   VectorXd a1 = AA->row(2*sample_idx);    
  //   VectorXd a2 = AA->row(2*sample_idx+1);
  //   VectorXd ci = cc->col(sample_idx);
  //   double b1 = (*bb)(0, sample_idx);
  //   double b2 = (*bb)(1, sample_idx);
  //   double di = (*dd)(sample_idx);    
  //   VectorXd x1 =  xData->row(4*sample_idx);
  //   VectorXd x2 =  xData->row(4*sample_idx+1);
  //   VectorXd x3 =  xData->row(4*sample_idx+2);
  //   VectorXd x4 =  xData->row(4*sample_idx+4);
  //   double y1 = (*yData)[4*sample_idx];
  //   double y2 = (*yData)[4*sample_idx+1];
  //   double y3 = (*yData)[4*sample_idx+2];
  //   double y4 = (*yData)[4*sample_idx+3];

  //   cout << a1.transpose()*theta + b1 << "\n";
  //   cout << a2.transpose()*theta + b2 << "\n";
  //   cout << ci.transpose()*theta + di<< "\n";
	
    
  //   bool isinls = QuasiconvexTools::isInlier(*xData, *yData, theta, sample_idx);
    
  // }
 
  
}

void HomogEstimator::store_model(int modelIndex, int numInliers)
{
  for (int i = 0; i < 9; ++i)
    final_model_params_[i] = *(models_denorm_[modelIndex]+i);
}

} // namespace usac

#endif
