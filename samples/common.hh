#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H
#include <string>
#include <libconfig.h++>
#include "usac.hh"
#include "problems/FundEstimator.hh"

namespace usac {

class ConfigFileReader: public libconfig::Config
{
 public:
  template <class T>
  bool getTopLevelParameterValue(const std::string settingName, const std::string parameterName,
                                 T& value) const
  {
    const libconfig::Setting& root = getRoot();
    try
    {
      const libconfig::Setting &setting = root[settingName.c_str()];
      if (!setting.lookupValue(parameterName.c_str(), value))
      {
        std::cout << settingName << ": parameter " << parameterName << " not found" << std::endl;
        return false;
      }
    } // end try
    catch(...)
    {
      std::cerr << settingName << " block not found; recheck configuration file" << std::endl;
      return false;
    }
    return true;
  }

  template <class T>
  bool getSubLevelParameterValue(const std::string settingName, const std::string subSettingName,
                                 const std::string parameterName, T& value) const
  {
    const libconfig::Setting& root = getRoot();
    try
    {
      const libconfig::Setting &setting = root[settingName.c_str()][subSettingName.c_str()];
      if (!setting.lookupValue(parameterName.c_str(), value))
      {
        std::cout << settingName << ":" << subSettingName
                  << ": parameter " << parameterName << " not found" << std::endl;
        return false;
      }
    }
    catch(...)
    {
      std::cerr << settingName << ":" << subSettingName << " block not found; recheck configuration file" << std::endl;
      return false;
    }
    return true;
  }

};

bool readParamsFromConfigFile(const ConfigFileReader& cfr, USACParams *p, double *inlier_th);

bool readParamsFromConfigFile(const ConfigFileReader& cfr, USACParams *cfg, FundParams *fund, double *inlier_th);

void print_usac_results(const USACResults &ret);

} // namespace usac

#endif
