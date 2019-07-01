#include "tractography.h"
#include "utilities.h"

#include <itkMacro.h> // needed for ITK_VERSION_MAJOR
#if ITK_VERSION_MAJOR < 5
#include "itkMultiThreader.h"
#else
#include "itkMultiThreaderBase.h"
#endif

#include "UKFTractographyCLP.h"

// TODO make configurable?
static const bool verbose = true;

namespace {
void ukf_setAndTell(ukfPrecisionType & x, const ukfPrecisionType y, const std::string & name)
{
  if (verbose) {
    x = y;
    std::cout << "- " << name << ": " << y << std::endl;
  }
}

void ukf_tell(const ukfPrecisionType & x, const std::string &name)
{
  if (verbose) {
    std::cout << "* " << name << ": " << x << std::endl;
  }
}
}; // anonymous namespace

int ukf_parse_cli(int argc, char** argv, UKFSettings& s)
{
  PARSE_ARGS ;

  /* Begin deprecation section */
  {
  /*
  *  Check for and print warning about invalid parameter `minGA`
  *
  *  GA is no longer used as a tracking threshold.
  *
  *  Please see the following for more information:
  *   https://github.com/pnlbwh/ukftractography/pull/64
  *   https://github.com/pnlbwh/ukftractography/pull/75
  *
  */

    // check infeasible default value as work-around because
    // Slicer CLI generates command line arguments for all
    // parameters, and the CLP doesn't indicate no-arg

    if (minGAArg.isSet() && minGA != 10000)
      {
      std::cerr << "Error: the `minGA` parameter is no longer valid because GA is not used! Please use 'stoppingThreshold' instead! Please see `--help` for more information." << std::endl;
      return EXIT_FAILURE;
      }
  }
  /* End deprecation section */

  ukfPrecisionType l_stoppingFA = stoppingFA;
  ukfPrecisionType l_stoppingThreshold = stoppingThreshold;
  ukfPrecisionType l_stepLength = stepLength;
  ukfPrecisionType l_recordLength = recordLength;
  ukfPrecisionType l_maxHalfFiberLength = maxHalfFiberLength;
  ukfPrecisionType l_seedingThreshold = seedingThreshold;
  ukfPrecisionType l_Qm = Qm;
  ukfPrecisionType l_Ql = Ql;
  ukfPrecisionType l_Qw = Qw;
  ukfPrecisionType l_Qkappa = Qkappa;
  ukfPrecisionType l_Qvic = Qvic;
  ukfPrecisionType l_Rs = Rs;
  ukfPrecisionType l_maxBranchingAngle = maxBranchingAngle;
  ukfPrecisionType l_minBranchingAngle = minBranchingAngle;

  // If sigmaSignal is not set minimum of voxel size is used for interpolation
  ukfPrecisionType SIGMA_SIGNAL = sigmaSignal;

  // HANDLE ERRORNOUS INPUT
  if (dwiFile.empty() || maskFile.empty() || tracts.empty()) {
    std::cout << "Error! Must indicate DWI data, mask and tracts output files!" << std::endl << std::endl ;
    return 1 ;	//This is to indicate that the module returns with error
  }

  if (numTensor == 1) {
    tractsWithSecondTensor.clear() ;	//Reassure the string is empty
  }

  if (l_maxHalfFiberLength <= 0) {
    std::cout << "Invalid maximum half fiber length!" << std::endl ;
    return 1 ;
  }

  if (std::ceil(l_maxHalfFiberLength / l_stepLength) <= 1) {
    std::cout << "Too large step length or too small fiber cutoff limit!" << std::endl ;
    return 1 ;
  }

  if (!freeWater && recordFreeWater) {
    std::cout << "In Order to use the \"--recordFreeWater\" flag, also the \"--freeWater\" flag must be used!" << std::endl;
    return 1 ;
  }

  if (freeWater && numTensor == 3) {
    std::cout << "In Order to use the free water estimation the number of Tensors (\"--numTensor\") must be set to 1 or 2.";
    std::cout << "(3-Tensor case not implemented yet.)" << std::endl;
    return 1 ;
  }

  if (l_recordLength < l_stepLength) {
    std::cout << "recordLength should be greater than stepLength" << std::endl;
    return 1 ;
  }

  // SETTING THE DEFAULT PARAMETERS
  std::string strModel = fullTensorModel ? "full model" : "simple model";
  std::string strFreeWater = freeWater ? " with free water estimation" : "";

  std::cout << "Using the " << numTensor << "T " << strModel << strFreeWater << ". Setting the default parameters accordingly:\n";
  std::cout << "\"*\": set by user\n";
  std::cout << "\"-\": default setting\n";

  if (seedsFile.empty()) {
    l_maxBranchingAngle = 0.0;
  }

  if (labels.size() == 0) {
    labels.push_back(1) ;	//Default to use label 1
  }

  if (l_stoppingFA == 0.15) {
    ukf_setAndTell(l_stoppingFA, l_stoppingFA, "stoppingFA");
  } else {
    ukf_tell(l_stoppingFA, "stoppingFA");
  }

  if (l_seedingThreshold == 0.0) {
    ukf_setAndTell(l_seedingThreshold, FULL_BRAIN_MEAN_SIGNAL_MIN,
               "seedingThreshold");  // Used to default to 2 times the FA threshold.
  } else {
    ukf_tell(l_seedingThreshold, "seedingThreshold");
  }

  if(!noddi){
    if(recordVic || recordKappa || recordViso){
      std::cout << "Can use recordVic or recordKappa or recordViso parameters only with noddi model";
      return EXIT_FAILURE;
    }
  }

  if (l_Qm == 0.0) {
    if (noddi){
      if (numTensor == 1)
        ukf_setAndTell(l_Qm, 0.0025, "Qm");
      else
        ukf_setAndTell(l_Qm, 0.001, "Qm");
    } else if (numTensor == 1) {
        ukf_setAndTell(l_Qm, 0.005, "Qm");//l_Qm = 0.0015;
    } else {
      if (fullTensorModel) {
        ukf_setAndTell(l_Qm, 0.002, "Qm");//l_Qm = 0.002;
      } else {
        ukf_setAndTell(l_Qm, 0.001, "Qm");//l_Qm = 0.001; was 0.003, changed to 0.001 for new Interp3Signal
      }
    }
  } else {
    ukf_tell(l_Qm, "Qm");
  }


  if (noddi){
    if ( l_Qkappa == 0.0)
      ukf_setAndTell (l_Qkappa, 0.01, "Qkappa");
    else
      ukf_tell(l_Qkappa, "Qkappa");
  }
  else {
    if (l_Ql == 0.0) {
      if (numTensor == 1) {
        ukf_setAndTell(l_Ql, 300.0, "Ql");//l_Ql = 25.0;
      } else if (numTensor == 2) {
        ukf_setAndTell(l_Ql, 50.0, "Ql");//was l_Ql = 100.0; for old Interp3Signal
      } else if (numTensor == 3) {
        ukf_setAndTell(l_Ql, 100.0, "Ql");//l_Ql = 150.0;
      }
    } else {
        ukf_tell(l_Ql, "Ql");
    }
  }


  if (l_Rs == 0.0) {
    if (numTensor == 1) {
      ukf_setAndTell(l_Rs, 0.01, "Rs");//l_Rs = 0.02;
    } else {
      if (fullTensorModel) {
        ukf_setAndTell(l_Rs, 0.01, "Rs");// = 0.01;
      } else {
        ukf_setAndTell(l_Rs, 0.02, "Rs");//was l_Rs = 0.015;for old Interp3Signal
      }
    }
  } else {
    ukf_tell(l_Rs, "Rs");
  }

  if (l_stepLength == 0.0) {
    if (numTensor == 1) {
      ukf_setAndTell(l_stepLength, 0.3, "stepLength");
    } else if (numTensor == 2) {
      ukf_setAndTell(l_stepLength, 0.3, "stepLength"); //was 0.2 for old Interp3Signal
    } else { // 3T
      ukf_setAndTell(l_stepLength, 0.15, "stepLength");
    }
  } else {
    ukf_tell(l_stepLength, "stepLength");
  }

  if (l_recordLength == 0.0) {
    if (numTensor == 1) {
      ukf_setAndTell(l_recordLength, 0.9, "recordLength");
    } else if (numTensor == 2) {
      ukf_setAndTell(l_recordLength, 0.9, "recordLength"); //was 0.2 for old Interp3Signal
    } else { // 3T
      ukf_setAndTell(l_recordLength, 0.45, "recordLength");
    }
  } else {
    ukf_tell(l_recordLength, "recordLength");
  }

  if (noddi){
    if (l_Qvic == 0.0)
      if (numTensor == 1)
        ukf_setAndTell(l_Qvic, 0.0005, "Qvic = Qviso");
      else
        ukf_setAndTell(l_Qvic, 0.004, "Qvic = Qviso");
    else
      ukf_tell(l_Qvic, "Qvic = Qviso");
  }
  else
    if (freeWater) {
      if (l_Qw == 0.0) {
        if (numTensor == 1) {
          ukf_setAndTell(l_Qw, 0.0025, "Qw"); // estimated in a paramsearch // 0.0025
        } else if (numTensor == 2) {
          ukf_setAndTell(l_Qw, 0.0015, "Qw"); // 0.0015
        }
      } else {
          ukf_tell(l_Qw, "Qw");
      }
    }

  ukf_tell(l_stoppingThreshold, "stoppingThreshold");

  if (seedsPerVoxel == 1) {
    std::cout << "- seedsPerVoxel: " << seedsPerVoxel << std::endl;
  } else {
    std::cout << "* seedsPerVoxel: " << seedsPerVoxel << std::endl;
  }

  // initializing settings
  //UKFSettings& s -- from function argument
    {
    s.record_fa = recordFA;
    s.record_nmse = recordNMSE;
    s.record_trace = recordTrace;
    s.record_state = recordState;
    s.record_cov = recordCovariance;
    s.record_free_water = recordFreeWater;
    s.record_tensors = recordTensors;
    s.record_Vic = recordVic;
    s.record_kappa = recordKappa;
    s.record_Viso = recordViso;
    s.transform_position = true; // TODO hard-coded :/
    s.store_glyphs = storeGlyphs;
    s.branches_only = false; // TODO hard-coded :/
    s.fa_min = l_stoppingFA;
    s.mean_signal_min = l_stoppingThreshold;
    s.seeding_threshold = l_seedingThreshold;
    s.num_tensors = numTensor;;
    s.seeds_per_voxel = seedsPerVoxel;
    s.min_branching_angle = l_minBranchingAngle;
    s.max_branching_angle = l_maxBranchingAngle;
    s.is_full_model = fullTensorModel;
    s.free_water = freeWater;
    s.noddi = noddi;
    s.stepLength = l_stepLength;
    s.recordLength = l_recordLength;
    s.maxHalfFiberLength = maxHalfFiberLength;
    s.labels = labels;
    s.num_threads = numThreads;

    s.Qm = l_Qm;
    s.Ql = l_Ql;
    s.Qw = l_Qw;
    s.Qkappa = l_Qkappa;
    s.Qvic = l_Qvic;
    s.Rs = l_Rs;

    // TODO these should be header-initialized once we use C++11
    s.p0 = P0;
    s.sigma_signal = SIGMA_SIGNAL;
    s.sigma_mask = SIGMA_MASK;
    s.min_radius = MIN_RADIUS;

    s.output_file = tracts;
    s.output_file_with_second_tensor = tractsWithSecondTensor;
    s.dwiFile = dwiFile;
    s.seedsFile = seedsFile;
    s.maskFile = maskFile;
    s.writeAsciiTracts = writeAsciiTracts;
    s.writeUncompressedTracts = writeUncompressedTracts;
    }

  return EXIT_SUCCESS;
}
