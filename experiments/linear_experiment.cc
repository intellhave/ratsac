#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
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
#include "extra.hh"
#include "LinearFit.hh"
#include "quasiconvexFunctions.hh"

using namespace std;
using namespace usac;
using namespace LinearTools;
using namespace MCTS;
using namespace ToolBox;

DEFINE_int32(N,100, "Data size");
DEFINE_int32(d,8, "Data dimensions");
DEFINE_double(th, 0.06, "Inlier Threshold");
DEFINE_double(p, 0.2, "Outlier Percentage");
DEFINE_double(insig, 0.1, "Inlier Variance");
DEFINE_double(osig, 1, "Outlier Variance");


string datasetFolder = "dataset/";
string dataConfigFile = "cfg.txt";
string resultFolder = "results/";

string datasetName = "data";
string configFile = "linear.cfg";

int N = FLAGS_N;
int d = FLAGS_d;
double insig = FLAGS_insig;
double osig = FLAGS_osig;
double inlier_th = FLAGS_th;
double outlierP = FLAGS_p;
int outlierStep = 5;
int outlierMin = 0;
int outlierMax = 30;
int nRansacRun = 10;


void genDataset(const Config &cfg){

  Eigen::MatrixXd xData;
  Eigen::VectorXd yData;
  Eigen::VectorXd thetaD;

  string datasetConfigFileName = datasetFolder + datasetName + "_config.txt";

  FILE *fconf = fopen(datasetConfigFileName.c_str(), "w");
  fprintf(fconf, "%d %d %lf %d %d %d",cfg.gen_N, cfg.gen_d, cfg.th, cfg.minOutlierRate, cfg.maxOutlierRate, cfg.outlierStep);
  fclose(fconf);
   
  for (int oP = cfg.minOutlierRate; oP < cfg.maxOutlierRate; oP+=cfg.outlierStep){

    string datasetFile = datasetName +  "_" + std::to_string(oP) + ".txt";

    srand(time(0));
    genRandomData(cfg.gen_N, cfg.gen_d, cfg.gen_insig, cfg.gen_osig, oP, xData, yData, thetaD);
    
    //write Generated Data to File
    ofstream dataFile(datasetFolder + datasetFile);
    for (int i=0; i<xData.rows(); i++){
      for (int j=0; j<xData.cols(); j++){
	dataFile << xData(i,j) <<" ";
      }
      dataFile << yData[i];
      dataFile << endl;
    }
    dataFile.close();
  }
}


void readDataset(int N, int d,int oP, Eigen::MatrixXd &xData, Eigen::VectorXd &yData){

  string datasetFile = datasetFolder + datasetName + "_" + std::to_string(oP) +".txt";
  FILE *finp;
  finp = fopen(datasetFile.c_str(),"r");

  xData.resize(N,d);
  yData.resize(N);

  for (int i=0; i<N; i++){
    double dt; 
    for (int j=0; j<d;j++){
      fscanf(finp, "%lf", &dt);
      xData(i,j) = dt;
    }
    fscanf(finp, "%lf", &dt) ;
    yData(i) = dt;
  }
  fclose(finp);
}

void runExp(Config &cfg){
  
  std::vector<string>methodList;
  string datasetConfigFileName = cfg.datasetFolder + cfg.datasetName + "_config.txt";

  /*Read Files*/
  FILE *fconf = fopen(datasetConfigFileName.c_str(), "r"); // Read dataset config file
  fscanf(fconf, "%d %d %lf %d %d %d", &N, &d, &inlier_th, &outlierMin, &outlierMax, &outlierStep);
  char method[20];
  
  /*Extract list of methods*/
  std::set<string> methods;
  int count = 0;
  while (!feof(fconf)){
    fscanf(fconf, "%s", method);
    string strMethod = method;
    if (strMethod.length()<2)
      continue;
    methods.insert(strMethod);
    if (methods.size()>count){
      count++;
      methodList.push_back(strMethod);
    }    
  }
  
  fclose(fconf); //Close the dataset config file

  methodList.clear();
  if (methodList.empty()){
    methodList = cfg.methodList;
  }
  
  /*Start*/
  for (int op = outlierMin; op < outlierMax; op+=outlierStep){
    cout<<"================================================\n";
    cout << op << "\n";
    //Read data from dataset
    Eigen::MatrixXd xData;
    Eigen::VectorXd yData;
    readDataset(N, d, op, xData, yData);    
    std::vector<int> methodInls;
    std::vector<double> methodRuntime;
   
    for (auto it=methodList.begin(); it!=methodList.end(); it++){
     
      string mt = *it;
      string rstFileName = resultFolder + datasetName +"_" + std::to_string(op)+ "_" +  mt +  ".txt";
      ofstream rstFile(rstFileName);
      cout  << "OutputFileName = " << rstFileName <<"\n";
      int nRun = 1;
      if (mt.compare("RANSAC")==0 || mt.compare("LORANSAC")==0 || mt.compare("MLESAC")==0 || mt.compare("PROSAC")==0 || mt.compare("MCTS")==0)
	nRun = cfg.nRansacRun;
      
      std::vector<int> inlsArr;
      std::vector<double> runTimeArr;      
      for (int run=0; run<nRun ; run++){
	//Initialize Linear Fitter
	LinearFit linearFitter(&xData, &yData, &cfg);
	linearFitter.sampleFile = "samples.txt";    /*File for storing samples picked during RANSAC*/
	//Read score for PROSAC
	std::vector<int> sorted;
	if (cfg.problem==FUNDMATRIX){
	  QuasiconvexTools::readSortedScoresFromFile(cfg.pathToScoreFile, sorted);
	  linearFitter.sorted = &sorted;
	}   
	//Start searching
	linearFitter.Fit(*it); /*Execute linear fit for the selected method*/
	
 	cout << "Run " <<run<< " Inls = " << linearFitter.inls << " Run Time = " << linearFitter.runTime << "\n";;
	//Write to file      
	rstFile << linearFitter.inls << "\t" << linearFitter.runTime <<"\n";
	
	inlsArr.push_back(linearFitter.inls);
	runTimeArr.push_back(linearFitter.runTime);
	
      }      
      rstFile.close();

      int tInls = round(std::accumulate(inlsArr.begin(), inlsArr.end(), 0.0)/nRun);
      double tRunTime = std::accumulate(runTimeArr.begin(), runTimeArr.end(), 0.0)/nRun;

      methodInls.push_back(tInls);
      methodRuntime.push_back(tRunTime);  
    }

    /*Write to Summary File*/
    string summFileName = cfg.resultFolder + cfg.datasetName + "_" + std::to_string(op) + "_all.txt";
    ofstream summFile(summFileName);
    for (int i=0; i<methodList.size(); i++){
      std::cout << op << "\t"<< methodList[i] << "\t" << methodInls[i] << "\t" << methodRuntime[i] <<"\n";
      summFile <<  methodList[i] << "\t" << methodInls[i] << "\t" << methodRuntime[i] <<"\n";       
    }
    summFile.close();
  }  
}

int main (int argc, char **argv){
  //  testLinearProgramming();
  
  if (argc < 2){
    cout << "Wrong parameters \n";
    return 0;
  }
          
  string configFileName = argv[2];
  Config cfg(configFileName);  
  cfg.readLinearGeneratorConfig();
  
  string type = argv[1];
  
  datasetFolder = cfg.datasetFolder;
  datasetName = cfg.datasetName;

  
  if (type.compare("-g")==0){
    std::cout << "=======GENERATING DATASET=============\n";
    genDataset(cfg);
  }
  else if (type.compare("-r")==0){
    runExp(cfg);
  }

}
