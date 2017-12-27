/*Test Homography Fitting */
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <numeric>

#include <gflags/gflags.h>
#include <cstdlib>
#include <ctime>
#include <set>
#include <Eigen/Core>

#include "problems/LinearEstimator.hh"
#include "common.hh"
#include "usac.hh"
#include "Chebyshev.hh"
#include "LinearFunctions.hh"
#include "config.hh"
#include "quasiconvexFunctions.hh"
#include "QuasiconvexFit.hh"
using namespace std;
using namespace usac;
//using namespace MCTS;
using namespace ToolBox;


int main(int argc, char **argv){

  string configFileName = argv[1];
  //string configFileName = "/home/intellhave/Dropbox/usac/build/homoexp.cfg";

  Config cfg(configFileName);

  Eigen::MatrixXd u1, u2; //Data

  QuasiconvexTools::readImageFeaturePointsFromFile(cfg.pathToDataFile, u1, u2);
  std::vector<int> sorted;
  QuasiconvexTools::readSortedScoresFromFile(cfg.pathToScoreFile, sorted);

  QuasiconvexFit hmFit(u1, u2, &sorted, &cfg);

  std::vector<int> methodInls;
  std::vector<double> methodRuntime;
  
  for (auto it = cfg.methodList.begin(); it!=cfg.methodList.end(); it++){
    string mt = *it;
    std::cout <<"Running  "<< mt << "\n";
    
    string rstFileName = cfg.resultFolder + cfg.datasetName + "_0_" + mt + ".txt";
    ofstream rstFile(rstFileName);
    int nRun = 1;

    if (mt.compare("RANSAC")==0 || mt.compare("LORANSAC")==0 || mt.compare("MLESAC")==0 || mt.compare("PROSAC")==0 || mt.compare("MCTS")==0)
      nRun = cfg.nRansacRun;

    std::vector<int> inlsArr;
    std::vector<double> runTimeArr;
    for (int run=0; run<nRun; run++){
      
      hmFit.Fit(*it);
      cout << "Run "  << run  << " Inls = " <<hmFit.inls << " Run Time = " << hmFit.runTime << "\n";
      //Write to file
      rstFile << hmFit.inls << "\t" << hmFit.runTime << "\n";
      
      inlsArr.push_back(hmFit.inls);
      runTimeArr.push_back(hmFit.runTime);    
    }
    rstFile.close();
    
    int tInls = round(std::accumulate(inlsArr.begin(), inlsArr.end(), 0.0)/nRun);
    double tRunTime = std::accumulate(runTimeArr.begin(), runTimeArr.end(), 0.0)/nRun;
    
    methodInls.push_back(tInls);
    methodRuntime.push_back(tRunTime);
  }

  /*Write to summary File*/
  string summFileName = cfg.resultFolder + cfg.datasetName + "_all.txt";
  ofstream summFile(summFileName);
  
  
  for (int i=0; i<cfg.methodList.size(); i++){
    std::cout << cfg.methodList[i] << "\t" << methodInls[i] << "\t" << methodRuntime[i] <<"\n";
    summFile << cfg.methodList[i] << "\t" << methodInls[i] << "\t" << methodRuntime[i] <<"\n";    
  }

  summFile.close();
  return 0;
}
