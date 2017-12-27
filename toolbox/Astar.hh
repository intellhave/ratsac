#ifndef ASTAR_HH
#define ASTAR_HH

#include <vector>
#include <set>
#include <Eigen/Core>
#include <queue>
#include <unordered_set>
#include <math.h>
#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "State.hh"
#include "extra.hh"


using namespace MCTS;
using namespace std;
using namespace Eigen;
using namespace LinearTools;
using namespace ToolBox;

class TNode{ //Tree node
public:
  Eigen::VectorXd p;
  std::set<int> b;
  double w;
  int f;
  int h;
  bool solved;
  std::set<int> v;

  bool operator < (const TNode &a) const  {
    return f > a.f;
  }
};


class Astar{
private:
  Eigen::MatrixXd *xData;
  Eigen::VectorXd *yData;
  double th;
  int N;
  int d; 

  priority_queue<TNode> Q;
  unordered_set<std::string> vHash;  
  int xnum;  

  std::set<int> oneToN = Tools::setOfInteger(0, N);
  std::set<int> empty; 
  

  int heuristic(set<int>H_in, set<int> chsub_old,  Eigen::VectorXd x00, //inp
		double &w, Eigen::VectorXd &x0_f, set<int> &chsub_old_f, set<int>& sub_f);  
  
  set<int> compute_upper(const std::set<int>&H, const Eigen::VectorXd &theta);
public:
  
 Astar(Eigen::MatrixXd *xData_, Eigen::VectorXd *yData_, double th_);
  
  
  Eigen::VectorXd pk;
  double wk;
  set<int> vk;
  set<int> bk;
  
  void startSearch();
  
  
  
};



#endif
