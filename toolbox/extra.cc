#include "extra.hh"
#include <iostream>
using namespace std;
using namespace Eigen;

namespace ToolBox{
  
  std::set<int> Tools::setOfInteger(int startIdx, int endIdx){
    set<int> S;
    for (int i=startIdx; i<endIdx; i++)
      S.insert(i);

    return S;
  }

  std::set<int> Tools::removeAtIndexes(const set<int>& s, const set<int>& indexes){

    std::set<int> r;
    vector<bool> flag;
    
    for (int i=0; i<s.size(); i++)
      flag.push_back(true);

    for (auto it=indexes.begin(); it!=indexes.end(); it++){
      flag[*it] = false;      
    }

    r.clear();
    for (int i=0; i<s.size(); i++){
      if (flag[i]){
	auto it = s.begin();
	advance(it, i);
	r.insert(*it);
      }
    }

    return r;
  }

  std::set<int> Tools::getSubsetAtIndexes(const std::set<int>&s, const std::set<int>& indexes){

    std::set<int> r;
    for (auto it = indexes.begin(); it!=indexes.end(); it++){
      auto its = s.begin();
      advance(its, *it);
      r.insert(*its);
    }
    
    return r;
  }

  void Tools::getDataSubset(const Eigen::MatrixXd *xData, const Eigen::VectorXd *yData, const set<int>&indexes,
		     Eigen::MatrixXd &sX, Eigen::VectorXd &sY){
    
    sX.resize(indexes.size(), xData->cols());
    sY.resize(indexes.size());
    int count = 0;
    
    for (auto it=indexes.begin(); it!=indexes.end(); it++){      
      sX.row(count) = (*xData).row(*it);
      sY[count] = (*yData)[*it];
      count++;
    }
    
  }

  int Tools::findInliers(const MatrixXd &xData, const VectorXd &yData, int nc, const VectorXd &theta, double th){
      
    int nConstraints = xData.rows();    
    int inls = 0;
    for (int i=0; i<nConstraints; i+=nc){
      bool isInls = true;
      for (int j=i; j<i+nc; j++){
	double r = xData.row(j)*theta - yData(j);
	//	cout << r << "\t";
	if (r>th){
	  isInls = false;
	  break;
	}
      }
      //      cout << "\n";
      if (isInls) inls++;
    }    
    return inls;
  }
  
   
}
