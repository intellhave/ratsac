/*Class to perform MCTS Search*/

#ifndef MSEARCH_HH
#define MSEARCH_HH

#include <cstdlib>
#include <Eigen/Core>
#include <vector>
#include <set>
#include <ctime>

#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "State.hh"
#include "tree.hh"
#include "mcts.hh"



namespace MCTS{
  class MSearch{
  private:

    //Linear Constraints
    Eigen::MatrixXd *xData;
    Eigen::VectorXd *yData;    
    int N; /*Number of constraints*/
    int nPoints; /*Number of points*/
    int dim; /*Problem dimensions*/
    double th;
    pType problemType;
    int treeTravelIter;

  public:
    /*Matrices for Quaisconvex*/
    Eigen::MatrixXd *AA;
    Eigen::MatrixXd *bb;
    Eigen::MatrixXd *cc;
    Eigen::VectorXd *dd;

    MSearch(Eigen::MatrixXd *x_, Eigen::VectorXd *y_, const int &nP_, const double &th_);
    
    Eigen::VectorXd theta; /*Solution*/
    double minimax; /*Minimax value of the solution*/
    int inls;
    std::set<int> basics; /*Basic of the solution*/
    double runTime;
    void startSearch();
       
  };
}





#endif
