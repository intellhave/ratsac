#ifndef  EXTRA_H
#define  EXTRA_H

#include <vector>
#include <set>
#include <string>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;
namespace ToolBox{

  class Tools{
  public:
    
    static inline std::vector<string> string_split(std::string str, char deliminator){
      std::vector<string> splitted;
      splitted.clear();
      stringstream ss(str);
      string se;
      
      while (getline(ss, se, deliminator)){
	splitted.push_back(se);
      }
      
      return splitted;     
    };
    
    
    static  std::set<int> setOfInteger(int startIdx, int endIdx);
    static  std::set<int> removeAtIndexes(const std::set<int>& s, const std::set<int>& indexes);
    static  std::set<int> getSubsetAtIndexes(const std::set<int>& s, const std::set<int>&indexes);
    
    static   void getDataSubset(const Eigen::MatrixXd *xData, const Eigen::VectorXd *yData, const set<int>& indexes,
			      Eigen::MatrixXd & sX, Eigen::VectorXd &sY);
    
    static   int findInliers(const MatrixXd &xData, const VectorXd &yData, int nc, const VectorXd &theta, double th);
  };
  
  

  
    
  }
    

#endif
