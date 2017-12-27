
#include "config.hh"

using namespace std;

namespace ToolBox{

  Config::Config(const string &cfgFile){
    configFile = cfgFile;
    cfgReader.readFile(configFile.c_str());
    readConfigFromFile();
    ransacConfigFile = cfgFile;
  }

  void Config::readConfigFromFile(){
    cfgReader.lookupValue("dataset_folder", datasetFolder);
    cfgReader.lookupValue("dataset_name", datasetName);
    cfgReader.lookupValue("result_folder", resultFolder);
    cfgReader.lookupValue("ransac_run", nRansacRun);
    cfgReader.lookupValue("min_samples", minRansacSamples);
    
    pathToDataFile = datasetFolder + datasetName + ".txt";
    pathToScoreFile = datasetFolder + datasetName + "_sorted.txt";    
    cfgReader.lookupValue("inlier_th", th);
    std::string methods;
    cfgReader.lookupValue("methods", methods);
    methodList = Tools::string_split(methods, ';');
    string strProblem;
    cfgReader.lookupValue("problem", strProblem);
    if (strProblem.compare("homography")==0)
      problem = HOMOGRAPHY;
    else if (strProblem.compare("affine")==0)
      problem = AFFINE;
    else if (strProblem.compare("fundmatrix")==0){
      problem = FUNDMATRIX;
      fmatrix = true;
    }
    else{
      problem = LINEAR;
      fmatrix = false;
    }
  }

  void Config::readLinearGeneratorConfig(){
    cfgReader.lookupValue("gen_insig", gen_insig);
    cfgReader.lookupValue("gen_osig", gen_osig);
    cfgReader.lookupValue("min_outlier_rate", minOutlierRate);
    cfgReader.lookupValue("max_outlier_rate", maxOutlierRate);
    cfgReader.lookupValue("outlier_step", outlierStep);
    cfgReader.lookupValue("N", gen_N);
    cfgReader.lookupValue("d", gen_d);
    //cfgReader.lookupValue("fmatrix", fmatrix);
    N = gen_N;
    d = gen_d;
  }


}
