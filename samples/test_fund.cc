#include <glog/logging.h>
#include <gflags/gflags.h>
#include "problems/FundEstimator.hh"
#include "common.hh"
#include "usac.hh"
#include <armadillo>


using namespace std;
using namespace usac;
using namespace arma;

DEFINE_bool(use_prosac, true, "");
DEFINE_string(p, "", "");

int main(int argc, char **argv) {
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  FundEstimator sac;
  ConfigFileReader cfr;
  cfr.readFile(argv[1]);
  FundParams fund_params;
  USACParams usac_params;
  double inlier_th, inl_th2;
  CHECK(readParamsFromConfigFile(cfr, &usac_params, &fund_params, &inlier_th));
  inl_th2 = inlier_th * inlier_th;

  string data_path = string(argv[2]);

  vector<int> sorting;
  ifstream sort_ifs(data_path + "/sorting.txt");
  CHECK(sort_ifs.is_open());
  string line;
  while (getline(sort_ifs, line)) {
    stringstream ss;
    ss << line;
    int idx;
    ss >> idx;
    sorting.push_back(idx);
  }
  sort_ifs.close();

  ifstream data_ifs(data_path + "/orig_pts.txt");
  CHECK(data_ifs.is_open());
  getline(data_ifs, line);
  stringstream ss;
  ss << line;
  int nr;
  ss >> nr;
  vector<double> pts;
  for (int i = 0; i < nr; i++) {
    double u0, v0, u1, v1;
    data_ifs >> u0 >> v0 >> u1 >> v1;
    pts.push_back(u0);
    pts.push_back(v0);
    pts.push_back(1);
    pts.push_back(u1);
    pts.push_back(v1);
    pts.push_back(1);
  }
  sac.init_params(usac_params);
  if (FLAGS_use_prosac)
    sac.init_data(nr, inl_th2, &sorting);
  else
    sac.init_data(nr, inl_th2);

  cerr << pts[0] << " " << pts[1] << " " << pts[2] << endl;
  sac.init_problem(fund_params, nr, &pts[0]);

  sac.solve();
  print_usac_results(sac.results_);

  ofstream ofs(data_path + "/"+FLAGS_p+"F.txt");
  // colmajor
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      ofs << sac.final_model_params_[i*3+j] << " ";
    ofs << std::endl;
  }

  ofstream inlier_ofs(data_path + "/"+FLAGS_p+"inliers.txt");
  for (auto i: sac.results_.inlier_flags_) {
    if (i)
      inlier_ofs << 1 << " ";
    else
      inlier_ofs << 0 << " ";
  }
}
