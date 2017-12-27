#ifndef CHEBYSHEV
#define CHEBYSHEV

#include <stdio.h>
#include <stdlib.h>
//#include <armadillo>
#include <set>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/SVD>

using namespace std;
//using namespace arma;


namespace Chebyshev{
  /*vec Descend(const mat &x, const mat &y, double &th, const vec &theta0,
    set<int>& solutionBasics, double &minimax);*/

  
  Eigen::VectorXd Descend(const Eigen::MatrixXd &xData, const Eigen::VectorXd &yData, double &th, const Eigen::VectorXd theta0, set<int>& solutionBasics, double &minimax);

  Eigen::VectorXd Descend(const Eigen::MatrixXd &xData, const Eigen::VectorXd &yData, double &th, const Eigen::VectorXd theta0, set<int>& solutionBasics, double &minimax, bool abs);

  
  Eigen::MatrixXd piv(const Eigen::MatrixXd &mat);
    
    
    
}



#endif
