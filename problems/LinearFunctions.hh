#ifndef LINEARTOOLS_HH
#define LINEARTOOLS_HH

#include <cstdlib>
//#include <armadillo>
#include <Eigen/Core>


//using namespace arma;

namespace LinearTools{

  //int findInliers(const mat &x, const vec &y, const vec &theta, double th);
  
  int getInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th,
		 Eigen::MatrixXd &xNew, Eigen::VectorXd &yNew);
  
  int findInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th);
  int findInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th, bool abs);
  
  
  //void genRandomData(const int N, const int d, const double insig, const double osig, const double outlierP, mat &x, vec &y, vec &theta);
  void genRandomData(const int N, const int d, const double insig, const double osig, const double outlierP, Eigen::MatrixXd &xData, Eigen::VectorXd &yData, Eigen::VectorXd &theta);
}



#endif
