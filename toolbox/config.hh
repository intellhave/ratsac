#ifndef CONFIG_HH
#define CONFIG_HH

#include <stdio.h>
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <libconfig.h++>
#include <vector>
#include "extra.hh"

using namespace std;

namespace ToolBox{

  enum prob{HOMOGRAPHY, AFFINE, FUNDMATRIX, LINEAR};
  
  class Config{
  private:
    libconfig::Config cfgReader;
    
  public:
    string configFile;
    string datasetFolder;
    string datasetName;
    string resultFolder;

    int nRansacRun;
    std::vector<string> methodList;
    
    double th;
    string pathToDataFile;
    string pathToScoreFile;
    string ransacConfigFile;
    bool fmatrix;

    prob problem;

    
    double gen_insig, gen_osig;
    int minOutlierRate, maxOutlierRate, outlierStep;
    int minRansacSamples;
    int gen_N, gen_d;
    int N, d; 
    Config(const string &cfgFile);
    void readConfigFromFile();
    void readLinearGeneratorConfig();
   };  
}
  

#endif
