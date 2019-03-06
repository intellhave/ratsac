#include "usac.hh"
#include "common.hh"

namespace usac {

bool readParamsFromConfigFile(const ConfigFileReader& cfr, USACParams *p, double *inlier_th)
{
  // read in parameters
  try
  {
    // first get all common parameters
    std::string rand_sampling_method, verif_method, local_opt_method, est_problem;
    if( !( cfr.getTopLevelParameterValue("common", "ransac_conf", p->conf_th) &&
           cfr.getTopLevelParameterValue("common", "min_sample_size", p->min_sample_size) &&
           cfr.getTopLevelParameterValue("common", "inlier_threshold", *inlier_th) &&
           cfr.getTopLevelParameterValue("common", "max_hypotheses", p->max_hypotheses) &&
           cfr.getTopLevelParameterValue("common", "max_solutions_per_sample", p->max_solns_per_sample) &&
           cfr.getTopLevelParameterValue("common", "prevalidate_sample", p->need_prevalidate_sample) &&
           cfr.getTopLevelParameterValue("common", "prevalidate_model", p->need_prevalidate_model) &&
           cfr.getTopLevelParameterValue("common", "degen_check", p->need_test_degeneracy) &&
           cfr.getTopLevelParameterValue("common", "random_sampling_method", rand_sampling_method) &&
           cfr.getTopLevelParameterValue("common", "verif_method", verif_method) &&
           cfr.getTopLevelParameterValue("common", "local_opt", local_opt_method)
	   )
	) //if
    {
      return false;
    }
    else
    {
      // verify parameter values
      if (p->conf_th < 0 || p->conf_th > 1)
      {
        std::cout << "RANSAC confidence value must be between 0 and 1" << std::endl;
        return false;
      }

      if (rand_sampling_method.compare("UNIFORM") == 0)
        p->sampling_method = SAMP_UNIFORM;
      else if (rand_sampling_method.compare("PROSAC") == 0)
        p->sampling_method = SAMP_PROSAC;
      else if (rand_sampling_method.compare("FILE")==0)
	p->sampling_method = SAMP_FILE;
	       
      else
      {
        std::cerr << "Random sampling method " << rand_sampling_method << " not recognized" << std::endl;
        return false;
      }

      if (verif_method.compare("STANDARD") == 0)
        p->verif_method = VERIF_STANDARD;
      else if (verif_method.compare("SPRT") == 0)
        p->verif_method = VERIF_SPRT;
      else
      {
        std::cerr << "Verification method " << verif_method << " not recognized" << std::endl;
        return false;
      }

      if (local_opt_method.compare("NO_LO") == 0)
        p->need_local_optim = false;
      else if (local_opt_method.compare("LOSAC") == 0)
        p->need_local_optim = true;
      else
      {
        std::cerr << "Local optimization method " << local_opt_method << " not recognized" << std::endl;
        return false;
      }
    }

    // read in PROSAC parameters if required
    if (p->sampling_method == SAMP_PROSAC)
    {
      if (!( cfr.getTopLevelParameterValue("prosac", "max_prosac_samples", p->prosac_max_samples) &&
             cfr.getTopLevelParameterValue("prosac", "beta", p->prosac_beta) &&
             cfr.getTopLevelParameterValue("prosac", "non_rand_conf", p->prosac_non_rand_conf) &&
             cfr.getTopLevelParameterValue("prosac", "min_stopping_length", p->prosac_min_stop_len)))
        // && cfr.getTopLevelParameterValue("prosac", "sorted_points_path", p->sortedPointsFile)))
      {
        return false;
      }
    }

    // read in SPRT parameters if required
    if (p->verif_method == VERIF_SPRT) {
      if ( !(cfr.getTopLevelParameterValue("sprt", "time_model", p->sprt_tM) &&
             cfr.getTopLevelParameterValue("sprt", "models_per_sample", p->sprt_mS) &&
             cfr.getTopLevelParameterValue("sprt", "delta", p->sprt_delta) &&
             cfr.getTopLevelParameterValue("sprt", "epsilon", p->sprt_epsilon)) )
      {
        return false;
      }
    }

    // read in local optimization parameters if required
    if (p->need_local_optim)
    {
      if ( !(cfr.getTopLevelParameterValue("losac", "inner_sample_size", p->losac_inner_sample_size) &&
             cfr.getTopLevelParameterValue("losac", "inner_ransac_repetitions", p->losac_inner_ransac_reps) &&
             cfr.getTopLevelParameterValue("losac", "threshold_multiplier", p->losac_threshold_multiplier) &&
             cfr.getTopLevelParameterValue("losac", "num_steps", p->losac_num_iterative_steps)) )
      {
        return false;
      }
    }
  }
  catch(...)
  {
    return false;
  }

  return true;
}

bool readParamsFromConfigFile(const ConfigFileReader& cfr, USACParams *p, FundParams *fund, double *inlier_th)
{
  // now read in parameters
  try
  {
    // usac parameters
    if (!readParamsFromConfigFile(cfr, p, inlier_th))
    {

      std::cerr << "Error reading USAC parameters from config file." << std::endl;
      return false;
    }

    if (p->need_test_degeneracy )
    {
      if ( !(cfr.getTopLevelParameterValue("problem_specific", "h_degen_threshold", fund->h_degen_th) &&
             cfr.getTopLevelParameterValue("problem_specific", "max_upgrade_samples", fund->max_upgrade_samples) ))
      {
        return false;
      }
    }

    std::string decomp_type;
    if ( !cfr.getTopLevelParameterValue("problem_specific", "matrix_decomposition", decomp_type) )
    {
      return false;
    }
    if (decomp_type.compare("LU") == 0)
    {
      fund->decomp_alg = DECOMP_LU;
    }
    else if (decomp_type.compare("QR") == 0)
    {
      fund->decomp_alg = DECOMP_QR;
    }
    else
    {
      std::cerr << "Matrix decomposition " << decomp_type << " not supported" << std::endl;
      return false;
    }
  }
  catch(...)
  {
    return false;
  }

  return true;
}

void print_usac_results(const USACResults &ret)
{
  std::cout << "Number of hypotheses/models: " << ret.hyp_nr_ << "/" << ret.model_nr_ << std::endl;
  std::cout << "Number of samples rejected by pre-validation: " << ret.rejected_sample_nr_ << std::endl;
  std::cout << "Number of models rejected by pre-validation: " << ret.rejected_model_nr_ << std::endl;
  std::cout << "Number of verifications per model: " <<
      (double) ret.total_points_verified_/(ret.model_nr_-ret.rejected_model_nr_) << std::endl;
  std::cout << "Max inliers/total points: " << ret.best_inlier_nr_ << "/" << ret.data_size << std::endl;
}

} //namespace usac
