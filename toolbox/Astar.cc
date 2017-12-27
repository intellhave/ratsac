
#include "Astar.hh"

using namespace std;
using namespace Eigen;
using namespace Chebyshev;
using namespace ToolBox;

Astar::Astar(Eigen::MatrixXd *xData_, Eigen::VectorXd *yData_,  double th_){
  xData = xData_;
  yData = yData_;
  
  th = th_;
  vHash.clear();
  N = (*xData).rows();
  d = (*xData).cols();
}


set<int> Astar::compute_upper(const set<int>&H, const Eigen::VectorXd &theta){

  set<int> s; s.clear();
  for (auto it = H.begin(); it!=H.end(); it++){
    double r = (*xData).row(*it)*theta - (*yData)[*it];
    if (fabs(r)>th)
      s.insert(*it);
  }    
  
  return s;
}

void Astar::startSearch(){
  xnum = 1;
  TNode node;
  
  Eigen::VectorXd theta0 = Eigen::VectorXd::Random(d); //Random init for Chebyshev

  double w;
  Eigen::VectorXd p;
  std::set<int> b;  
  std::set<int> bs;
  oneToN = Tools::setOfInteger(0, N);
  empty.clear();

  std::set<int> temp1;
  std::set<int> temp2;

  p = Chebyshev::Descend(*xData, *yData, th, theta0, b, w);
  node.p = p;
  node.b = b;
  node.w = w;
  node.v.clear();
  
  set<int> v = compute_upper(oneToN,node.p);
  int UN_cnt = v.size();
  
  node.h = heuristic(oneToN, empty, node.p, w, p, temp1, temp2);
  node.f = node.v.size() + node.h;
  
  if (node.w <=th){
    pk = node.p;
    wk = node.w;
    vk = node.v;
    bk = node.b;    
    return;
  }
  node.solved = false;   
  int minOUT = UN_cnt;
  Q.push(node);

  //  priority_queue<TNode> tempQueue;
  //tempQueue.push(node);    
  
  
  while (true){

    node = Q.top();
    //cout << "+ " << node.v.size() << "-----" <<std::flush;
    if (node.solved){
      pk = node.p;
      wk = node.w;
      vk = node.v;
      bk = node.b;    
      return;
    }

    //    cout << node.f << "\t" << node.h <<"\n";
    
    Q.pop();
    
    for (auto it=node.b.begin(); it!=node.b.end(); it++){
      TNode child;
      child.v = node.v;
      child.v.insert(*it);
      
      string vstr = MCTS::stringFromSet(child.v);
      unordered_set<string>::iterator fit = vHash.find(vstr);
      if (fit!=vHash.end()){
	continue;
      }
      
      set<int> H = Tools::setOfInteger(0,N); 
      H = Tools::removeAtIndexes(H, child.v); //H is the support set of the current node

      Eigen::MatrixXd xH;  Eigen::VectorXd yH;
      Tools::getDataSubset(xData, yData, H, xH, yH);
      std::set<int> bs;
      
      child.p = Chebyshev::Descend(xH, yH, th, node.p, bs, child.w);
      child.b = Tools::getSubsetAtIndexes(H, bs);
      //Compute heuristics;
      child.h = heuristic(H, empty, child.p, w, p, temp1, temp2);
      child.h = max(child.h, node.h-1);

      child.h = 0;
      
      child.f = child.v.size()+child.h;

      vHash.insert(vstr);

      v = compute_upper(H, p);
      UN_cnt = v.size();
      
      if (minOUT > UN_cnt + child.v.size())
	minOUT = UN_cnt + child.v.size();

      if (child.h == UN_cnt && child.f <=minOUT){ //h is lower bound, UN_cnt is upper bound. 
	child.solved = true;
	child.p = p;
	child.w = w;	
	set<int> Hv = Tools::getSubsetAtIndexes(H, v);
	for (auto ithv = Hv.begin(); ithv!=Hv.end();  ithv++)
	  child.v.insert(*ithv);	  
      }
      else{
	child.solved = 0;
	if (child.w <=th){
	  pk = child.p;
	  wk = child.w;
	  vk = child.v;	  
	  return;
	}
      }     
      Q.push(child);      
    }           
  }              
}

int Astar::heuristic(set<int>H_in, set<int> chsub_old, Eigen::VectorXd x00,
		     double &w, Eigen::VectorXd &x0_f, set<int> &chsub_old_f, set<int>& sub_f){


  //Generate data for the support set
  Eigen::MatrixXd xH; Eigen::VectorXd yH;
  Tools::getDataSubset(xData, yData, H_in, xH, yH);

  //Perform Chebyshev fitting
  double val;
  set<int> chsub;
  Eigen::VectorXd x0 = Chebyshev::Descend(xH, yH, th, x00, chsub, val);
  
  //TEST
  //Eigen::VectorXd res = xH*x0 - yH;
  //cout << res <<"\n";
  //cout << "Val = " << val << "\t";
  //for (auto itt = chsub.begin(); itt!=chsub.end(); itt++)
  //  cout << res(*itt) << "\t" << std::flush;

  set<int> H_chsub = Tools::getSubsetAtIndexes(H_in, chsub);
  set<int> H = Tools::removeAtIndexes(H_in, chsub); /*Remove all basics for the current set*/
  
  /*Check the removed points, see if it can be added back to the support set*/
  for (auto it = chsub_old.begin(); it!=chsub_old.end(); it++){
    double r = (*xData).row(*it)*x0 - (*yData)[*it]; 
    if (fabs(r)<val)
      H.insert(*it);          
    else
      chsub_old_f.insert(*it);         
  }
  
  int cyc;  /*Heuristic value*/
  set<int> chsub_new;    
  if (val>th){       /*If the current fit is not feasible*/
    set<int> H1 = H;
    cyc = heuristic(H1, H_chsub, x0, w, x0_f, chsub_new, sub_f);    
  }
  else {
    x0_f = x0;
    w = val;
    sub_f = H;
    sub_f.insert(H_chsub.begin(), H_chsub.end());
    return 0;
  }
  
  /*Iterate through removed points*/
  for (auto it=chsub_new.begin(); it!=chsub_new.end(); it++){
    sub_f.insert(*it);
    if (sub_f.size()<=chsub.size()){
      cyc++;
      return cyc;;
    }

    double r = (*xData).row(*it)*x0_f - (*yData)[*it];
    if (fabs(r) < th)
      continue;

    double valn;
    set<int> bs;
    Eigen::MatrixXd xF; Eigen::VectorXd yF;
    /*Add back the removed point, perform fitting for the new set*/
    Tools::getDataSubset(xData, yData, sub_f, xF, yF);   
    Eigen::VectorXd x1 = Chebyshev::Descend(xF, yF, th, x0_f, bs, valn);

    if (valn > th){
      sub_f = Tools::removeAtIndexes(sub_f, bs);
      if (sub_f.size() <= chsub.size()){
	cyc++;
	return cyc;
      }
      cyc++;
    }
    else{
      x0_f = x1;
      w = valn;
    }
  }      
}
  



//DEBUG - check 
    // if (tempQueue.size() > 50){
    //   cout << "----------------\n";
    //   while (tempQueue.size()>0){
    // 	TNode tempNode = tempQueue.top();
    // 	std::cout << tempNode.f << "\t" <<std::flush;
    // 	tempQueue.pop();
    //   }
    // }
    
