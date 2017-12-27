#ifndef QUASICONVEX_FUNCTIONS
#define QUASICONVEX_FUNCTIONS

#include <stdio.h>
#include <cstdlib>
#include <Eigen/Core>
#include <vector>
#include <set>
#include "common.hh"
#include <ctime>


using namespace std;
using namespace Eigen;

namespace QuasiconvexTools{

  extern void genLinearMatrixFromQuasiconvex(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, double th, MatrixXd &xData, VectorXd &yData);

  extern void readImageFeaturePointsFromFile(std::string fileName, MatrixXd  &u1, MatrixXd &u2);

  extern void readSortedScoresFromFile(std::string fileName, std::vector<int>&sorted);

  extern bool isInlier(const MatrixXd &xData, const VectorXd &yData, const VectorXd &theta, int ptsIdx);

  extern double computeError(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, int ptsIdx);

  extern  int findInliers(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, double th);

  extern int findInliers(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, double th, std::set<int> &inlierSet, std::set<int>&outlierSet, std::set<int>& basicSet);

  extern int findInliersFromLinearConstraints(const MatrixXd &xData, const VectorXd &yData, const VectorXd &theta);



}


#endif
