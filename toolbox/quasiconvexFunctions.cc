#include "quasiconvexFunctions.hh"

namespace QuasiconvexTools{

  /****************************************************************************************/
  void readImageFeaturePointsFromFile(std::string fileName, Eigen::MatrixXd  &u1, Eigen::MatrixXd  &u2){

    int nPoints;
    FILE *finp = fopen(fileName.c_str(), "r");
    fscanf(finp, "%d", &nPoints);
    u1.resize(2, nPoints);
    u2.resize(2, nPoints);
    int pt = 0;

    while (!feof(finp) && pt < nPoints){
      double ux, uy, vx, vy;
      fscanf(finp, "%lf %lf %lf %lf", &ux, &uy, &vx, &vy);
      u1(0,pt) = ux; u1(1,pt) = uy;
      u2(0,pt) = vx; u2(1,pt) = vy;
      pt++;
    }

    fclose(finp);
  }

  /****************************************************************************************/
  void readSortedScoresFromFile(std::string fileName, std::vector<int>&sorted){
    int nPoints;
    FILE *finp = fopen(fileName.c_str(), "r");

    fscanf(finp, "%d", &nPoints);
    sorted.clear();
    for (int i=0; i<nPoints; i++){
      int pts;
      fscanf(finp, "%d", &pts);
      sorted.push_back(pts);
    }
    fclose(finp);
  }


  /****************************************************************************************/
  void genLinearMatrixFromQuasiconvex(const MatrixXd &A, const MatrixXd &b, const MatrixXd  &c, const VectorXd &d, double th, MatrixXd &xData, VectorXd &yData){

    int N = A.rows();
    xData.resize(2*N, 8);
    yData.resize(2*N);

    int rowIdx = 0;
    for (int i=0; i<N-1; i=i+2){
      VectorXd a1 = A.row(i);
      VectorXd a2 = A.row(i+1);
      int idx = round(i/2);

      double b1 = b(0, idx);
      double b2 = b(1, idx);

      double thd = th*d(idx);
      VectorXd thc = th*c.col(idx);

      xData.row(rowIdx) = a1 + a2 - thc; yData(rowIdx) = thd - b1 - b2; rowIdx++;
      xData.row(rowIdx) = a2 - a1 - thc; yData(rowIdx) = thd - b2 + b1; rowIdx++;
      xData.row(rowIdx) = a1 - a2 - thc; yData(rowIdx) = thd - b1 + b2; rowIdx++;
      xData.row(rowIdx) = -a1 - a2- thc; yData(rowIdx) = thd + b1 + b2; rowIdx++;
    }
  }

  /****************************************************************************************/
  bool isInlier(const MatrixXd &xData, const VectorXd &yData, const VectorXd &theta,  int ptsIdx){
    int nc = 4;/*nc : number of constraints for each data point*/
    int startIdx = ptsIdx * nc;
    bool inls = true;
    for (int i=0; i<nc; i++){
      double r = xData.row(startIdx+i)*theta - yData[startIdx+i];
      cout << xData.row(startIdx+i) << "\n";
      cout << r <<"\n" << std::flush;
      if (r>0){
	return false;
      }
    }

    return true;

  }
  /****************************************************************************************/
  double computeError(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, int ptsIdx){

    VectorXd Axpb(2);
    Axpb(0) = A.row(2*ptsIdx)*theta + b(0, ptsIdx);
    Axpb(1) = A.row(2*ptsIdx+1)*theta + b(1, ptsIdx);
    double cxpd = c.col(ptsIdx).dot(theta) + d(ptsIdx);
    double r = Axpb.lpNorm<1>();
    r = r/cxpd;
    return fabs(r);

  }


  /****************************************************************************************/
  int findInliers(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, double th){

    int inls = 0;
    int nPoints = A.rows()/2;

    for (int i=0; i<nPoints; i++){
      double r = computeError(A, b, c, d, theta, i);
      if (fabs(r) <= th) inls++;
    }

    return inls;
  }

  /****************************************************************************************/
  int findInliers(const MatrixXd &A, const MatrixXd &b, const MatrixXd &c, const VectorXd &d, const VectorXd &theta, double th, std::set<int> &inlierSet, std::set<int>&outlierSet, std::set<int>& basicSet){

    int inls = 0;
    int nPoints = A.rows()/2;
    inlierSet.clear(); outlierSet.clear(); basicSet.clear();
    std::vector<double> residuals;
    for (int i=0; i<nPoints; i++){
      double r = computeError(A, b, c, d, theta, i);
      residuals.push_back(r);      
      if (fabs(r) <= th){
	inls++;
	inlierSet.insert(i);	
      }
      else {
	outlierSet.insert(i);
      }      
    }
    

    
    
  }
  
  /****************************************************************************************/

  int findInliersFromLinearConstraints(const MatrixXd &xData, const VectorXd &yData, const VectorXd &theta){

    int nPoints = (xData.rows()/4);
    int inls = 0;
    for (int i=0; i<nPoints; i++){
      if (isInlier(xData, yData, theta, i))
	inls++;
    }
    return inls;

  }


}
