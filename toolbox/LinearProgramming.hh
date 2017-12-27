#ifndef LINEAR_PROGRAMMING_HH
#define LINEAR_PROGRAMMING_HH

#include <stdio.h>
#include <cstdlib>
#include <Eigen/Core>

#include "OsiClpSolverInterface.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"



using namespace Eigen;

namespace ToolBox{

class LinearProgramming{

private:
  MatrixXd *A;
  VectorXd *c;
  VectorXd *b;
  VectorXd *var_lb;
  VectorXd *var_ub;
  int N, d;

public:
  LinearProgramming(MatrixXd *A_, VectorXd *b_, VectorXd *c_, VectorXd *var_lb_, VectorXd *var_ub_);
  void solve();
  VectorXd soln;
};

}



#endif
