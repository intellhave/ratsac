/*
  Perform Linear Fitting for a particular method: RANSAC, LORANSAC, MLESAC, PROSAC, L1, LINF, ASTAR
*/

#include <iostream>
#include "LinearFit.hh"

using namespace std;

namespace ToolBox{
  
  /*
  Constructor for the LinearFit class that solve
  max_sum_i( II(abs(X_i*theta - Y_i)<= th)) 
  where II(.) is the indicator function
  x_: pointer to an Eigen matrix that stores X data
  y: pointer to an Eigen vector that stores Y data
  cfg_: pointer to a Config Object that stores some parameters for the experiments
  
  */
  LinearFit::LinearFit(Eigen::MatrixXd *x_, Eigen::VectorXd *y_, Config *cfg_){
    
    xData = x_;
    yData = y_;
    cfg = cfg_;
    inlier_th = cfg->th;    
    N = xData->rows();
    d = xData->cols();
    bestTheta.resize(d);
    configFile = cfg_->configFile.c_str();
    sorted = NULL;
    
  }
  

  /*Perform Linearfit with the method selected*/
  void LinearFit::Fit(const std::string& method){
    
    cout << "--------------------------------------------------\n";
    cout << "Executing linear fit by using " << method << "\n";
    cout << "--------------------------------------------------\n";
    
    /*Set parameters for ransac variants*/
    sampling_method = SAMP_UNIFORM;
    verif_method = VERIF_STANDARD;
    need_local_optim = false;
    mle_sac = false;
    
    if (method.compare("RANSAC")==0){
      ransacFit();
      sac.writeSamplesToFile(sampleFile);
    }
    else if (method.compare("LORANSAC")==0){
      sampling_method = SAMP_FILE; /*Use the same sample set as RANSAC*/
      need_local_optim = true;     /*Set this flag in RASAC to enable LO-RANSAC*/      
      ransacFit();
    }
    else if(method.compare("MLESAC")==0){
      mle_sac = true;
      ransacFit();
    }
    else if (method.compare("PROSAC")==0){
      sampling_method = SAMP_PROSAC;
      ransacFit();
    }
    else if(method.compare("L1")==0){
      startTime = clock();
      L1Fit();
      stopTime = clock();
    }
    else if (method.compare("LINF")==0){
      startTime = clock();
      LInfFit();
      stopTime = clock();
      inls = LinearTools::findInliers(*xData, *yData, bestTheta, inlier_th);      
    }
    else if (method.compare("MCTS")==0){
      startTime = clock();      
      /*Initialization*/
      double kappa = 10.0;
      LInfFit();
      //inls = LinearTools::findInliers(*xData, *yData, bestTheta, inlier_th*3);
      Eigen::MatrixXd xNew; Eigen::VectorXd yNew;      
      int linfInls = LinearTools::getInliers(*xData, *yData, bestTheta, inlier_th*kappa, xNew, yNew);
      /*Point the data to the new matrix*/
      auto *xDataOld = xData; auto *yDataOld = yData;
      xData = &xNew;
      yData = &yNew;      
      N = xNew.rows();     
      MCTSFit();
      stopTime = clock();
      inls = LinearTools::findInliers(*xDataOld, *yDataOld, state::bestTheta, inlier_th);
      //int inls2 = LinearTools::findInliers(*xData, *yData, state::bestTheta, inlier_th);                 
    }
    else if (method.compare("ASTAR")==0){
      startTime = clock();
      Astar astarSearch = Astar(xData, yData, inlier_th);
      astarSearch.startSearch();      
      stopTime = clock();
      cout << N - astarSearch.vk.size() <<"\n";
      inls = LinearTools::findInliers(*xData, *yData, astarSearch.pk, inlier_th);           
    }        
    runTime  = (stopTime - startTime)*1.0/CLOCKS_PER_SEC;

    /*Find Inliers*/
    
  }
  


  void LinearFit::ransacFit(){
    //Read config file
    cfr.readFile(configFile);
    USACParams usac_params;
    LinearParams linear_params;
   
    double th1;
    CHECK(readParamsFromConfigFile(cfr, &usac_params, &th1));
    usac_params.min_hypotheses = cfg->minRansacSamples;
    usac_params.min_sample_size = d;
    usac_params.sampling_method = sampling_method;
    usac_params.need_local_optim = need_local_optim;
    usac_params.verif_method = verif_method;
    usac_params.mle_sac = mle_sac;
    sac.allSamples.clear();
    
    sac.init_params(usac_params);
    sac.init_data(N, inlier_th);
    sac.init_problem(linear_params, xData, yData);
        
    if (sampling_method == SAMP_FILE)
      sac.readSamplesFromFile(sampleFile);

    startTime = clock();
    sac.solve();
    stopTime = clock();
    print_usac_results(sac.results_);
    //inls = sac.results_.best_inlier_nr_;
    /*Recompute the inliers - Use the inliner computation for all methods*/
    inls = LinearTools::findInliers(*xData, *yData, sac.final_model_params_, inlier_th);
    //cout << sac.final_model_params_ << "\n";
  }

  
  void LinearFit::MCTSFit(){
    Eigen::VectorXd rand_theta = Eigen::VectorXd::Random(d);
    
    double minimax;
    Eigen::VectorXd thetaf = Chebyshev::Descend((*xData), (*yData), inlier_th, rand_theta, basics, minimax);
    
    set<int>vSet;  
    set<int>sSet;

    vSet.clear();
    for (int i=0; i<N; i++)
      sSet.insert(i);

    /*Static variables to store some common parameters during the Monte-Carlo Search process*/
    state::totalIter = 0; // Total number of iterations
    state::maxInls = 0;   // Maximum consensus size so far
    state::bestTheta.resize(d); //Best result so far

    if (state::stateMap.size()>0){
      for (auto it = state::stateMap.begin(); it!=state::stateMap.end(); it++){
	delete it->second;
      }
    }
    state::stateMap.clear(); //Hash table to stores fitted samples
    
    clock_t tick = clock();
    state initState(xData, yData, vSet, sSet, inlier_th, MCTS::DOWN, thetaf);
    //Create a tree with the init state;
    Tree *searchTree = new Tree(&initState);
    state::maxInls = 0;
    for (int i=0; i<1; i++){
      Node *currentNode = searchTree->root;      
      while (!currentNode->nodeState->terminated()){
	state *currentState = currentNode->nodeState;
	
	cout << currentState->minimax << " --minimax " << currentState->state_th << " --th "  <<endl;// << "  "<<currentNode<<endl;
	
	int maxChildren = currentState->direction==DOWN?
	  currentState->basics.size():currentState->violationSet.size();    
	int iter = maxChildren;
	iter  = 6;
	Node *prevNode = currentNode;

	currentNode = uctSearch(searchTree, currentNode, iter);

	if (prevNode == currentNode){
	  cout << "---------=======-------------"<<endl;
	  break;
	}
      }
    }
    //double mcts_time = (clock()-tick)*1.0/CLOCKS_PER_SEC;
    //inls = state::maxInls;
    bestTheta = state::bestTheta;
    //    inls = LinearTools::findInliers(*xData, *yData, state::bestTheta, inlier_th);
    //cout << state::bestTheta <<"\n";
  }


  void LinearFit::LInfFit(){    

    double maxRes = std::numeric_limits<double>::max();
    /*Initially include all constraints*/    
    std::set<int> constraintSet = Tools::setOfInteger(0, N);
    Eigen::MatrixXd xFit;
    Eigen::VectorXd yFit;

    Eigen::VectorXd theta = Eigen::VectorXd::Random(d);
    std::set<int> fBasics;
    
    while (maxRes > inlier_th){      
      Tools::getDataSubset(xData, yData, constraintSet, xFit, yFit);
      theta = Chebyshev::Descend(xFit, yFit, inlier_th, theta, fBasics, maxRes);      
      /*Remove points with largest residuals*/
      constraintSet = Tools::removeAtIndexes(constraintSet, fBasics);                  
    }
    bestTheta = theta;        
  } 

  void LinearFit::L1Fit(){
    Eigen::MatrixXd lA(2*N, d+N);
    Eigen::VectorXd lb(2*N);
    Eigen::VectorXd lc(d+N);
    Eigen::VectorXd var_lb(N+d);
    Eigen::VectorXd var_ub(N+d);
    /*Fill the matrix and vector for linear programming*/
    for (int i=0; i<N; i++){
      for (int j=0; j<d; j++) lA(i,j) = (*xData)(i,j);
      for (int j=0; j<N; j++) if (j==i) lA(i, j+d)=-1; else lA(i,j+d)=0;
      lb[i] = (*yData)[i] + inlier_th;      
    }
    
    for (int i=0; i<N; i++){
      for (int j=0; j<d; j++) lA(N+i, j) = -1.0*(*xData)(i,j);
      for (int j=0; j<N; j++) if (j==i) lA(N+i, d+j)=-1; else lA(N+i, d+j) = 0;
      lb[i+N] = -(*yData)[i] + inlier_th;      
    }
    

    for (int i=0; i<d; i++) {
      var_lb[i] = -std::numeric_limits<double>::max();
      var_ub[i] = std::numeric_limits<double>::max();
      lc[i] = 0;
    }

    for (int i=0; i<N; i++){
      var_lb[d+i] = 0;
      var_ub[d+i] = std::numeric_limits<double>::max();
      lc[d+i] = 1;
    }

    
    LinearProgramming lp_sover(&lA, &lb, &lc, &var_lb, &var_ub);
    lp_sover.solve();

    //cout << lp_sover.soln <<"\n";
    for (int i=0; i<d; i++){
      bestTheta[i] = lp_sover.soln[i];
    }
    
    inls = LinearTools::findInliers(*xData, *yData, bestTheta, inlier_th);
    
    
    
  }
  
  

  
}






