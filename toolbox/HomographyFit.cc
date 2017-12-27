#include "HomographyFit.hh"


namespace ToolBox{

  HomographyFit::HomographyFit(const Eigen::MatrixXd &u1_, const Eigen::MatrixXd &u2_, std::vector<int>* sorted_, Config *cfg_){
    u1 = u1_;
    u2 = u2_;
    cfg = cfg_;
    sorted = sorted_;

    inlier_th = cfg->th;
    nPoints = u1.cols();
    N = nPoints*2;
  }

  void HomographyFit::genMatrixHomography(){
    A.resize(N, 8);
    b.resize(2, nPoints);
    c.resize(8, nPoints);
    d.resize(nPoints);
    for (int i=0; i<nPoints; i++){

      double x1 = u1(0,i); double x2 = u1(1,i);
      double xp1 = u2(0,i); double xp2 = u2(1,i);

      Eigen::VectorXd ai1(8), ai2(8), ci(8);
      ai1 << x1, x2, 1, 0, 0, 0, -xp1*x1, -xp1*x2;
      ai2 << 0, 0, 0, x1, x2, 1, -xp2*x1, -xp2*x2;
      A.row(2*i) = ai1;
      A.row(2*i+1) = ai2;

      b(0,i) = -xp1;
      b(1,i) = -xp2;
      ci << 0,0,0,0,0,0,x1,x2;
      c.col(i) = ci;
      d(i)=1;
    }
  }

  void HomographyFit::ransacFit(){
    vector<double> pts;
    for (int i=0; i<nPoints; i++)
    {
      pts.push_back(u1(0,i)); pts.push_back(u1(1,i)); pts.push_back(1);
      pts.push_back(u2(0,i)); pts.push_back(u2(1,i)); pts.push_back(1);
    }

    HomogEstimator sac;
    double th_temp;
    USACParams usac_params;
    cout << cfg->ransacConfigFile <<"\n";
    cfr.readFile("homoexp.cfg");
    cfr.readFile(cfg->ransacConfigFile.c_str());
    CHECK(readParamsFromConfigFile(cfr, &usac_params, &th_temp));

    usac_params.min_hypotheses = cfg->minRansacSamples;
    cout << usac_params.min_hypotheses << "\n";
    usac_params.min_sample_size = 4;
    usac_params.sampling_method = sampling_method;
    usac_params.need_local_optim = need_local_optim;
    usac_params.verif_method = verif_method;
    usac_params.mle_sac = mle_sac;

    sac.allSamples.clear();
    sac.init_params(usac_params);
    sac.init_data(nPoints, inlier_th, sorted);
    sac.init_problem(nPoints, &pts[0]);

    genMatrixHomography();
    QuasiconvexTools::genLinearMatrixFromQuasiconvex(A, b, c, d, inlier_th, xData, yData);
    sac.AA = &A;
    sac.bb = &b;
    sac.cc = &c;
    sac.dd = &d;
    sac.xData = &xData;
    sac.yData = &yData;


    sac.solve();
    print_usac_results(sac.results_);

    //Test inlier counting functions
    Eigen::VectorXd ransacTheta(8);
    for (int i=0; i<8; i++)
      ransacTheta[i] = sac.final_model_params_[i]*1.0/sac.final_model_params_[8];
    inls = QuasiconvexTools::findInliers(A, b, c, d, ransacTheta, inlier_th);
    cout << inls << "\n";

  }
  
  void HomographyFit::LInfFit(){
    double maxRes = std::numeric_limits<double>::max();
    
    //std::set<int> pointSet = Tools::setOfInteger(0, nPoints);
    std::set<int> constraintSet = Tools::setOfInteger(0, N);
    std::set<int> fBasics;
    Eigen::VectorXd theta = Eigen::VectorXd::Random(8);    
    while (maxRes > 0.00000001){
      cout << maxRes << "\t" << constraintSet.size()<< "\n";
      // constraintSet.clear();
      // for (auto it = pointSet.begin(); it!=pointSet.end(); it++){
      // 	for (int rIdx=4*(*it); rIdx < 4*(*it)+4; rIdx++)
      // 	  constraintSet.insert(rIdx);
      // }
      Eigen::MatrixXd xFit;
      Eigen::VectorXd yFit;    
      Tools::getDataSubset(&xData, &yData, constraintSet, xFit, yFit);
      theta = Chebyshev::Descend(xFit, yFit, inlier_th, theta, fBasics, maxRes, false);
      std::set<int> basics; basics.insert(*fBasics.begin());
      constraintSet = Tools::removeAtIndexes(constraintSet, fBasics);           
    }
    
    bestTheta = theta;
    inls = QuasiconvexTools::findInliers(A, b, c, d, bestTheta, inlier_th);    
    
  }
  
  
  void HomographyFit::L1Fit(){
    
    Eigen::MatrixXd lA(xData.rows(), xData.cols()+nPoints);
    Eigen::VectorXd lb = yData;
    Eigen::VectorXd lc(xData.cols()+nPoints);
    Eigen::VectorXd var_lb(xData.cols() + nPoints);
    Eigen::VectorXd var_ub(xData.cols()+ nPoints);
    
    for (int i=0; i<xData.rows(); i++){
      for (int j=0; j<xData.cols(); j++) lA(i,j) = xData(i,j);
      for (int j=0; j<nPoints; j++)
	if (j==(i/4))
	  lA(i, xData.cols()+j)=-1;
	else
	  lA(i, xData.cols()+j)=0;      
    }

    for (int i=0; i<xData.cols(); i++){
      lc[i] = 0;
      var_lb[i] = -std::numeric_limits<double>::max();
      var_ub[i] = std::numeric_limits<double>::max();
    }
    for (int i=xData.cols(); i<xData.cols()+nPoints; i++){
      lc[i] = 1;
      var_lb[i] = 0;
      var_ub[i] = std::numeric_limits<double>::max();      
    }

    LinearProgramming lp_solve(&lA, &lb, &lc, &var_lb, &var_ub);
    lp_solve.solve();

    bestTheta.resize(xData.cols());
    for (int i=0; i<xData.cols(); i++)
      bestTheta[i] = lp_solve.soln[i];

    inls = QuasiconvexTools::findInliers(A, b, c, d, bestTheta, inlier_th);
  }
  
    
  
  void HomographyFit::Fit(const string &method)
  {
    startTime = clock();
    /*Prepare parameters for RANSAC variants*/
    sampling_method = SAMP_UNIFORM;
    verif_method  = VERIF_STANDARD;
    need_local_optim = false;
    mle_sac = false;

    if (method.compare("RANSAC")==0){
      ransacFit();
    }
    else if (method.compare("LORANSAC")==0){
      need_local_optim = true;
      ransacFit();
    }
    else if (method.compare("PROSAC")==0){
      sampling_method = SAMP_PROSAC;
      ransacFit();
    }
    else if (method.compare("MLESAC")==0){
      mle_sac = true;
      ransacFit();
    }
    else if (method.compare("LINF")==0){
      genMatrixHomography();
      QuasiconvexTools::genLinearMatrixFromQuasiconvex(A, b, c, d, inlier_th, xData, yData);              LInfFit();
    }
    else if (method.compare("L1")==0){
      genMatrixHomography();
      QuasiconvexTools::genLinearMatrixFromQuasiconvex(A, b, c, d, inlier_th, xData, yData);
      L1Fit();      
    }
    else if (method.compare("MCTS")==0){

      genMatrixHomography();
      QuasiconvexTools::genLinearMatrixFromQuasiconvex(A, b, c, d, inlier_th, xData, yData);
      MSearch homoMCTS(&xData, &yData, nPoints, inlier_th);
      /*Assign quasiconvex matrices to MSearch*/
      homoMCTS.AA = &A;
      homoMCTS.bb = &b;
      homoMCTS.cc = &c;
      homoMCTS.dd = &d;
      homoMCTS.startSearch();
      inls = QuasiconvexTools::findInliers(A, b, c, d, homoMCTS.theta, inlier_th);    

    }
    
    stopTime = clock();
    runTime = (stopTime - startTime)*1.0/CLOCKS_PER_SEC;


  }






}
