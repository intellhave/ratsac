#include <cstdlib>
#include <iostream>
#include <ctime>
#include <random>

#include "LinearFunctions.hh"

//using namespace arma;
namespace LinearTools{

  /*
  int findInliers(const mat &x, const vec &y, const vec &theta, double th){
    
    vec r = abs(x*theta - y);
    int inlsCount = 0;
    for (int i=0; i<r.n_elem; i++)
      if (r[i]<=th) inlsCount++;
    
    return inlsCount;
    
  }
  */

  /*Get inliers given theta and threshold*/
  int getInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th, 
		 Eigen::MatrixXd &xNew, Eigen::VectorXd &yNew)
  {
    /*Compute residuals*/
    Eigen::VectorXd r = (x*theta - y);
    //    if (abs)
    r = r.cwiseAbs();    
    int inlsCount = 0;
    int nPoints = x.rows();
    
    bool *isInls = new bool[r.size()];
    memset(isInls, 0, nPoints*sizeof(bool));    
	   
    for (int i=0; i<r.size(); i++){
      if (r[i]<=th)
	{
	  inlsCount++;
	  isInls[i] = 1;	  
	}
    }
    /*Now, start pushing inliers to the new matrix*/ 
    xNew.resize(inlsCount, x.cols());
    yNew.resize(inlsCount);
    int rowCount = 0;
    for (int i=0; i<r.size(); i++){ 
      if (isInls[i]){
	/*I'm busy to check so I did this stupid row copy*/
	rowCount++;
	for (int j=0; j<x.cols(); j++)
	  xNew(rowCount-1,j) = x(i,j);
	yNew[rowCount-1] = y[i];  	  
      }		  
    }       
    return inlsCount;
  }
  
  
	   
  int findInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th){
    return findInliers(x, y, theta, th, true);
  }
  
  
  int findInliers(const Eigen::MatrixXd &x, const Eigen::VectorXd &y, const Eigen::VectorXd &theta, double th, bool abs){
    
    Eigen::VectorXd r = (x*theta - y);
    if (abs)
      r = r.cwiseAbs();
    
    int inlsCount = 0;
    for (int i=0; i<r.size(); i++)
      if (r[i]<=th) inlsCount++;
    
    return inlsCount;
    
  }

  /*
  void genRandomData(const int N, const int d, const double insig, const double osig, const double outlierP, mat &x, vec &y, vec &theta){
    
    srand(time(0));
    arma_rng::set_seed(time(0));
    x = randn<mat>(N, d);  
    theta = randn<vec>(d);
    y = x*theta + insig*randn<vec>(N);

    int nOutlier;
    if (outlierP <1)
      nOutlier = round(N*outlierP);
    else
      nOutlier = round(N*outlierP*1.0/100.0);
	
    //Gen outlier
    for (int i=0; i<nOutlier; i++){
      double temp = osig*((double)randn()/RAND_MAX);
      if (temp>0) temp = temp+2*osig;
      y[i] = y[i] + temp;
    }             
  }
  */
  
  void genRandomData(const int N, const int d, const double insig, const double osig, const double outlierP, Eigen::MatrixXd &xData, Eigen::VectorXd &yData, Eigen::VectorXd &theta){

    /*Initialze random generators*/
    srand(time(0));
    std::default_random_engine generator;
    std::normal_distribution<double> inls_dist(0, insig);
    std::normal_distribution<double> out_dist(0, osig);
    
    //xData.setRandom(N,d);    
    //theta.setRandom(d);
    xData.resize(N,d);
    theta.resize(d);
    yData.resize(N);

    /*Generate random x data*/
    for (int i=0; i<N; i++){
      for (int j=0; j<d-1; j++){
	xData(i,j)=(double)rand()/RAND_MAX;
      }
      xData(i,d-1) = 1;
    }
    /*Generate random theta*/
    for (int i=0; i<d; i++)
      theta(i) = (double)rand()/RAND_MAX;

    //    cout << theta <<"\n";
    /*Add Gaussian noise*/
    for (int i=0; i<N; i++){
      //double gNoise = insig*((double)rand()/RAND_MAX);
      double gNoise = inls_dist(generator);
      yData[i] = xData.row(i)*theta + gNoise;      
    }
    //cout << yData << "\n";
    /*Add Outliers*/
    int nOutlier;
    if (outlierP <1)
      nOutlier = round(N*outlierP);
    else
      nOutlier = round(N*outlierP/100);
    
    for (int i=0; i<nOutlier; i++){
      //double outlier_noise = osig*((double)rand()/RAND_MAX);
      double outlier_noise = (out_dist(generator));

      double r = yData[i] - xData.row(i)*theta;
      // if (r > 0)
      // 	outlier_noise = 4*outlier_noise;
      // else
      // 	outlier_noise = -outlier_noise;
	
      yData[i] = yData[i] + outlier_noise;

      for (int j=0; j<d; j++){
      	//xData(i,j) += osig*(double)rand()/RAND_MAX;
	xData(i,j)+=out_dist(generator);
      }
      
    }             
  }   
}

