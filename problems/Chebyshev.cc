#define ARMA_DONT_PRINT_ERRORS

#include <cstdlib>
#include <iostream>
#include <armadillo>
#include "Chebyshev.hh"
#include <math.h>
#include <Eigen/QR>

using namespace std;
using namespace arma;
const double MAX_NUM = 10000000.0;

namespace Chebyshev{
  
  Eigen::VectorXd Descend(const Eigen::MatrixXd &xData, const Eigen::VectorXd &yData, double &th, const Eigen::VectorXd theta0, set<int>& solutionBasics, double &minimax){
    return Chebyshev::Descend(xData, yData, th, theta0, solutionBasics, minimax, true);
    
  }

  
  Eigen::VectorXd Descend(const Eigen::MatrixXd &xData, const Eigen::VectorXd &yData, double &th, const Eigen::VectorXd theta0, set<int>& solutionBasics, double &minimax, bool abs){

    Eigen::MatrixXd X;
    Eigen::VectorXd Y;
    
    if (!abs){
      X.resize(xData.rows(), xData.cols());
      Y.resize(yData.size());
      X << xData;
      Y << yData;
    }
    else {
      X.resize(2*xData.rows(), xData.cols());
      Y.resize(2*yData.size());
      X << xData, -xData;  
      Y << yData, -yData;      
    }


    int N = X.rows();       //Number of constraints
    int n = X.cols();       //Dimensions
    
    Eigen::VectorXd theta(n);
    if (!isnan(theta0(0)))
      theta << theta0;
    else
      theta = Eigen::VectorXd::Random(n);
    
    Eigen::VectorXd r = X*theta - Y;
    int maxIndexR;
    double d = r.maxCoeff(&maxIndexR);
    minimax = d;
    
    vector<int> basics;
    basics.clear();
    /* Find the inital set of basisc */
    for (int i=0; i<N; i++){
      if (fabs(d-r[i])<0.0000001){
	basics.push_back(i);
	r[i] =  d;
	if (basics.size()==n+1)
	  break;
      }
    }
    
    int basicSize = basics.size();
    Eigen::MatrixXd A(basicSize, n);
    
    // Extract baiscs to matrix A    
    for (int i=0; i<basicSize; i++){     
      A.row(i) =  X.row(basics[i]);
    }	
    
    // If number of points in the basics is less than dimension+1, continue to search    
    while (basicSize<n+1){

      Eigen::VectorXd t = Eigen::VectorXd::Ones(basicSize);
      Eigen::MatrixXd B = -1*A*A.transpose();
      //Eigen::VectorXd c = B.colPivHouseholderQr().solve(t);
      //Eigen::VectorXd c = B.partialPivLu().solve(t);
      Eigen::VectorXd c = B.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(t);
 
      //mat B = A*A.t();
      //vec c = solve(-B, t);
      //vec c = B\t;

      Eigen::VectorXd y = A.transpose()*c;
      Eigen::VectorXd dd = d*Eigen::VectorXd::Ones(N) - r;
      Eigen::VectorXd lambda(N);
      Eigen::VectorXd dd_den = X*y + Eigen::VectorXd::Ones(N);
      
      for (int i=0; i<N;i++){
	lambda[i] = dd[i]/dd_den[i];
	if (isnan(lambda[i]) || lambda[i]<=0) lambda[i] = MAX_NUM;
      }

      int j;
      double lambda_j = lambda.minCoeff(&j);      

      /*Add j into the basic set*/
      basicSize = basicSize + 1;      
      basics.push_back(j);      
      A.conservativeResize(basicSize, n);
      A.row(basicSize-1)= X.row(j);//      A = join_cols(A, X.row(j));
     
      theta = theta + lambda_j*y;
      r = (X*theta - Y);                                       // Recompute residuals
      d = r.maxCoeff(&maxIndexR);
      minimax = d;
      
      for (int i=0; i<basics.size(); i++)
      	r[basics[i]] = d;      
    }

    // std::cout << basics.size() << "\n";

    Eigen::MatrixXd CC(n+1,n+1);
    CC.col(0) = Eigen::VectorXd::Ones(n+1);
    for (int i=1; i<n+1; i++)
      CC.col(i) = A.col(i-1);
    
    //CC = CC.colPivHouseholderQr().inverse();
    //CC = CC.fullPivLu().inverse();
    CC = piv(CC);
	
    //CC = CC.inverse();
    Eigen::VectorXd C1 = CC.row(0).transpose();//    vec C1 = CC.row(0).t();    
    int totalSign = 0;    
    for (int i=0; i<C1.size(); i++)
      totalSign += C1[i]>=0?1:0;    
    
    int iter = 0;
    //    cout << endl;
    while (totalSign < n && iter<100){      
      iter++;
      int p1;
      double minval = C1.minCoeff(&p1);//      int p1 = C1.index_min();      
      Eigen::VectorXd y(CC.rows()-1);//    vec y(CC.n_rows-1);
      for (int i=0; i<CC.rows()-1; i++)
	y[i] = CC(i+1, p1)/CC(0, p1);
      
      Eigen::VectorXd dd = d*Eigen::VectorXd::Ones(N) - r;
      Eigen::VectorXd dd_den = X*y + Eigen::VectorXd::Ones(N);
      
      Eigen::VectorXd t(N);       //vec t = (d*ones<vec>(N)-r)/(X*y+ones<vec>(N));
      for (int i=0; i<N; i++){
	t[i] = dd[i]/dd_den[i];
	if (t[i]<=0 || isnan(t[i])) t[i] = MAX_NUM;
      }
      
      int j;
      double min_t = t.minCoeff(&j);//      int j  = t.index_min();
      
      theta = theta + min_t*y;
      
      Eigen::VectorXd t1(n+1);
      t1[0] = 1; for (int i=1; i<n+1; i++) t1[i] = X(j, i-1);      //      vec t1(n+1); t1[0] = 1; for (int i=1; i<n+1; i++) t1[i] = X(j, i-1);

      Eigen::MatrixXd lambda = t1.transpose()*CC;  //      mat lambda = (t1.t()*CC);     
      int beta = p1; 
      Eigen::VectorXd Cb = CC.col(beta)/lambda(beta); // vec Cb = CC.col(beta)/lambda(beta);      

      Eigen::MatrixXd L1 = Eigen::MatrixXd::Ones(n+1,1)*lambda;
      Eigen::MatrixXd L2 = Cb*Eigen::MatrixXd::Ones(1,n+1);
      CC = CC - L1.cwiseProduct(L2);   //      CC = CC - ( (ones<mat>(n+1,1)*lambda ) % ( Cb*ones<mat>(1,n+1)) ) ;
      CC.col(beta) = Cb;
      
      basics[beta] = j;
      
      r = (X*theta - Y);
      d = r.maxCoeff(&maxIndexR);//     d = max(r);
      minimax = d;
      
      for (int i=0; i<basics.size(); i++)
	r[basics[i]] = d;
    
      C1 = CC.row(0).transpose();
      totalSign = 0;
      for (int i=0; i<C1.size(); i++)
	totalSign += C1[i]>=0?1:0;          
    }
    
    int l = xData.rows();    
    solutionBasics.clear();
    
    for (int i=0; i<basics.size(); i++){
      if (basics[i]>=l)
	basics[i] -= l;
      solutionBasics.insert(basics[i]);
    }
    
    //    cout << "....descend"<<endl;
    return theta;    
  }

  

  Eigen::MatrixXd piv(const Eigen::MatrixXd &mat){

    double epsilon = 0.01;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon*std::max(mat.cols(), mat.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    
      // template<typename _Matrix_Type_>
    //   _Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
    // {
    //   Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    //   double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    //   return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
    // }
  }

  /*
  vec Descend(const mat &x_, const mat &y_, double &th, const vec &theta0,
	      set<int>&solutionBasics, double &minimax)

  {
    mat X = join_cols(x_, -x_);
    vec Y = join_cols(y_, -y_);

    int N = X.n_rows;       //Number of constraints
    int n = X.n_cols;       //Dimensions
    
    vec theta = theta0;     //Solution
    
    vec r = X*theta - Y;    
    double d = max(r);    
    
    vector<int> basics;
    basics.clear();
    // Find the inital set of basisc
    for (int i=0; i<N; i++){
      if (d-r[i] < 0.000001){
	basics.push_back(i);
	r[i] =  d;
      }
    }
    minimax = d;

    if (basics.size()==0)
      basics.push_back(index_max(r));
    
    int basicSize = basics.size();
    mat A(1,n);

    // Extract baiscs to matrix A    
    std::vector<int>::iterator it = basics.begin();
    A.row(0) = X.row(*it);
    for (int i=1; i<basicSize; i++){
      it++;
      if (it!=basics.end())	  
	A = join_cols(A, X.row(*it));
    }	

    // If number of points in the basics is less than dimension+1, continue to search
    while (basicSize<n+1){           
      vec t = ones<vec>(basicSize);
      mat B = A*A.t();
      vec c = solve(-B, t);
      //vec c = B\t;
      
      vec y = A.t()*c;
      vec dd = d*ones<vec>(N) - r;
      vec lambda = dd/(X*y+ ones<vec>(N));
      
      for (int i=0; i<N; i++){
	if (isnan(lambda[i])) lambda[i] = 0;
	if (lambda[i]<=0) lambda[i] = MAX_NUM;
      }
      
      double lambda_j = min(lambda);
      int j = lambda.index_min();
      
      basicSize = basicSize + 1;      
      basics.push_back(j);      
      A = join_cols(A, X.row(j));

      theta = theta + lambda_j*y;                 
      r = abs(X*theta - Y);                                       // Recompute residuals
      d = max(r);      
      
      for (int i=0; i<basics.size(); i++)
      	r[basics[i]] = d;      
    }       

    mat onevec = ones<mat>(A.n_rows,1);
    mat CC = join_rows(onevec, A);    
    CC = pinv(CC);
    vec C1 = CC.row(0).t();
    
    int totalSign = 0;
    
    for (int i=0; i<C1.n_elem; i++)
      totalSign += C1[i]>=0?1:0;    
    
    int iter = 0;
    while (totalSign < n && iter<100){
      iter++;
      double minval = min(C1);
      int p1 = C1.index_min();
      
      vec y(CC.n_rows-1);
      for (int i=0; i<CC.n_rows-1; i++)
	y[i] = CC(i+1, p1)/CC(0, p1);

      vec t = (d*ones<vec>(N)-r)/(X*y+ones<vec>(N));
      for (int i=0; i<N; i++){
	if (t[i]<=0) t[i] = MAX_NUM;
      }
      
      double min_t = min(t);
      int j  = t.index_min();

      //cout << CC <<endl;
      //cout<< CC.row(0) << endl;
      theta = theta + min_t*y;	    
      vec t1(n+1); t1[0] = 1; for (int i=1; i<n+1; i++) t1[i] = X(j, i-1);      
      mat lambda = (t1.t()*CC);     
      int beta = p1; 
      vec Cb = CC.col(beta)/lambda(beta);      
      
      CC = CC - ( (ones<mat>(n+1,1)*lambda ) % ( Cb*ones<mat>(1,n+1)) ) ;
      CC.col(beta) = Cb;
      
      basics[beta] = j;

      r = (X*theta - Y);
      d = max(r);
      minimax = d;
      
      for (int i=0; i<basics.size(); i++)
	r[basics[i]] = d;
    
      C1 = CC.row(0).t();
      totalSign = 0;
      for (int i=0; i<C1.n_elem; i++)
	if (C1[i]>=0)
	  totalSign++;	
    }
    
    int l = x_.n_rows;
    
    solutionBasics.clear();
    
    for (int i=0; i<basics.size(); i++){
      if (basics[i]>l)
	basics[i] -= l;
      solutionBasics.insert(basics[i]);
    }

    // for (std::vector<int>::iterator it = basics.begin(); it!=basics.end(); it++){
    //   if (*it > l)
    // 	*it = *it - l;
    //   solutionBasics.insert(*it);
    // }
    
    // for (int i=0; i<basics.size(); i++)
    //   if (basics[i]>l) basics[i] = basics[i]-l;
    return theta;
  }
  */
  
}





      // if (isnan(minval)){
      // 	cout << C1 <<endl;
      // 	cout << invMat<<endl;
      // 	cout << prevInv <<endl;
      // 	cout << A << endl;
      // }
