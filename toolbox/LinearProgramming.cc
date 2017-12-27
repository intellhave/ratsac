
#include "LinearProgramming.hh"

namespace ToolBox{

  LinearProgramming::LinearProgramming(MatrixXd *A_, VectorXd *b_, VectorXd *c_, VectorXd *var_lb_, VectorXd *var_ub_){
    A = A_; b = b_; c = c_; var_lb = var_lb_; var_ub = var_ub_;    
    N = A->rows();
    d = A->cols();    
    soln.resize(d);
    
    if (c->size() != d){
      std::cout << "Wrong inputs" << "\n";
    }    
  }
  
  void LinearProgramming::solve(){    
    OsiSolverInterface *si = new OsiClpSolverInterface();
    /*Copy data to Clp matrix and vectors*/
    int n_cols = d;
    double *objective = new double[n_cols];
    /*Copy Objective*/
    for (int i=0; i<d; i++){  objective[i] = (*c)[i];   }
    
    
    double *cols_lb = new double[n_cols];
    double *cols_ub = new double[n_cols];
    for (int i=0; i<d; i++){
      cols_lb[i] = (*var_lb)[i];
      cols_ub[i] = (*var_ub)[i];
    }
    

    int n_rows = N;
    double *rows_lb = new double[n_rows];
    double *rows_ub = new double[n_rows];
    
    CoinPackedMatrix *matrix = new CoinPackedMatrix(false, 0, 0);
    matrix->setDimensions(0, n_cols);
    
    for (int i=0; i<N; i++){
      CoinPackedVector row;
      for (int j=0; j<d; j++){ row.insert(j, (*A)(i,j)); }
      matrix->appendRow(row);
      rows_lb[i] = -1*si->getInfinity();
      rows_ub[i] = (*b)[i];      
    }
    
    si->loadProblem(*matrix, cols_lb, cols_ub, objective, rows_lb, rows_ub);
    si->initialSolve();
    if (si->isProvenOptimal()){
      
      const double *solution;
      solution = si->getColSolution();
      for (int i=0; i<d;i++)
	soln[i] = solution[i];
    }else  {
      std::cout <<"Failed to run LP\n";
      
    }
    
  }

  
}



