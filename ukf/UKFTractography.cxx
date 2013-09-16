/**
 * \file UKFTractography.cxx
 * \brief Main file of the project
 *
 * In this file the input arguments are parsed and passed on to the tractography
 * Also the model choice happens here
*/

#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "BRAINSThreadControl.h"
#include "filter_model.h"
#include "tractography.h"
#include "UKFTractographyCLP.h"


static const bool verbose = true;

void setAndTell(double & x, const double y, std::string name)
{
  if (verbose) {
    x = y;
    std::cout << "- " << name << ": " << y << std::endl;
  }
}

void tell(double & x, std::string name)
{
  if (verbose) {
    std::cout << "* " << name << ": " << x << std::endl;
  }
}

int main(int argc, char **argv)
{

  PARSE_ARGS ;

  std::cout << std::endl;

  // CONSTANTS
  bool FULL_BRAIN                 = false;
  const double SIGMA_SIGNAL       = 1.66;
  const double SIGMA_MASK 	      = 0.5;
  const double P0 			          = 0.01;
  const double MIN_RADIUS 		    = 0.87;
  const double FULL_BRAIN_GA_MIN	= 0.18;
  const double D_ISO			        = 0.003; // Diffusion coefficient of free water

  // NOTE:  When used as share libary one must be careful not to permanently reset number of threads
  //        for entire program (i.e. when used as a slicer modules.
  //        This also addresses the issue when the program is run as part of a batch processing
  //        system so that the number of cores allocated by scheduler is respected rather than
  //        blindly using all the cores that are found.
  //        This implementation is taken from extensive testing of the BRAINSTools
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numThreads);
  const int actuallNumThreadsUsed = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
  {
    std::cout << "Found " << actuallNumThreadsUsed << " cores on your system.\n";
    std::cout << "Running tractography with " << actuallNumThreadsUsed << " thread(s).\n";
  }

  std::cout << std::endl;

  // HANDLE ERRORNOUS INPUT
  if (dwiFile.empty() || maskFile.empty() || tracts.empty()) {
    std::cout << "Error! Must indicate DWI data, mask and tracts output files!" << std::endl << std::endl ;
    return 1 ;	//This is to indicate that the module returns with error
  }

  if (numTensor == 1) {
    tractsWithSecondTensor.clear() ;	//Reassure the string is empty
  }

  if (numTensor <= 0 || numTensor > 3) {
    std::cout << "Invalid tensor number!" << std::endl << std::endl ;
    return 1 ;
  }

  if (maxHalfFiberLength <= 0) {
    std::cout << "Invalid maximum half fiber length!" << std::endl ;
    return 1 ;
  }

//   if (stepLength <= 0){
//     std::cout << "Invalid step length!" << std::endl ;
//     return 1 ;
//   }

  if (std::ceil(maxHalfFiberLength / stepLength) <= 1) {
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


  // SETTING THE DEFAULT PARAMETERS
  std::string strModel = simpleTensorModel ? "simple model" : "full model";
  std::string strFreeWater = freeWater ? " with free water estimation" : "";

  std::cout << "Using the " << numTensor << "T " << strModel << strFreeWater << ". Setting the default parameters accordingly:\n";
  std::cout << "\"*\": set by user\n";
  std::cout << "\"-\": default setting\n";

  if (seedsFile.empty()) {
    FULL_BRAIN = true;
    maxBranchingAngle = 0.0;
  }

  if (labels.size() == 0) {
    labels.push_back(1) ;	//Default to use label 1
  }

  if (minFA == 0.15) {
    setAndTell(minFA, minFA, "minFA");
  } else {
    tell(minFA, "minFA");
  }

  if (seedFALimit == 0.0) {
    setAndTell(seedFALimit, minFA, "seedFALimit");  // Used to default to 2 times the FA threshold.
  } else {
    tell(seedFALimit, "seedFALimit");
  }

  if (Qm == 0.0) {
    if (numTensor == 1) {
      setAndTell(Qm, 0.005, "Qm");//Qm = 0.0015;
    } else {
      if (!simpleTensorModel) {
        setAndTell(Qm, 0.002, "Qm");//Qm = 0.002;
      } else {
        setAndTell(Qm, 0.003, "Qm");//Qm = 0.003;
      }
    }
  } else {
    tell(Qm, "Qm");
  }


  if (Ql == 0.0) {
    if (numTensor == 1) {
      setAndTell(Ql, 300.0, "Ql");//Ql = 25.0;
    } else if (numTensor == 2) {
      setAndTell(Ql, 100.0, "Ql");//Ql = 100.0;
    } else if (numTensor == 3) {
      setAndTell(Ql, 100.0, "Ql");//Ql = 150.0;
    }
  } else {
    tell(Ql, "Ql");
  }


  if (Rs == 0.0) {
    if (numTensor == 1) {
      setAndTell(Rs, 0.01, "Rs");//Rs = 0.02;
    } else {
      if (!simpleTensorModel) {
        setAndTell(Rs, 0.01, "Rs");// = 0.01;
      } else {
        setAndTell(Rs, 0.015, "Rs");//Rs = 0.015;
      }
    }
  } else {
    tell(Rs, "Rs");
  }

  if (stepLength == 0.0) {
    if (numTensor == 1) {
      setAndTell(stepLength, 0.3, "stepLength");
    } else if (numTensor == 2) {
      setAndTell(stepLength, 0.2, "stepLength");
    } else { // 3T
      setAndTell(stepLength, 0.15, "stepLength");
    }
  } else {
    tell(stepLength, "stepLength");
  }

  if (freeWater) {
    if (Qw == 0.0) {
      if (numTensor == 1) {
        setAndTell(Qw, 0.0025, "Qw"); // estimated in a paramsearch // 0.0025
      } else if (numTensor == 2) {
        setAndTell(Qw, 0.0015, "Qw"); // 0.0015
      }
    } else {
      tell(Qw, "Qw");
    }
  }

  tell(minGA, "minGA");

  if (seedsPerVoxel == 1) {
    std::cout << "- seedsPerVoxel: " << seedsPerVoxel << std::endl;
  } else {
    std::cout << "* seedsPerVoxel: " << seedsPerVoxel << std::endl;
  }

  if (normalizedDWIData) {
    outputNormalizedDWIData = false ;
  }

  if (weightsOnTensors.empty()) {
    for (int i = 0; i < numTensor; i++) {
      weightsOnTensors.push_back(1.0 / numTensor) ;
    }
  } else {
    if (static_cast<int>(weightsOnTensors.size()) != numTensor) {
      std::cout << "Wrong number of weights on tensors!" << std::endl << std::endl ;
      exit(1) ;
    }

    double weight_accumu = 0 ;
    for (int i = 0; i < numTensor; i++) {
      weight_accumu += weightsOnTensors[i] ;
    }
    if (std::abs(weight_accumu - 1.0) > 0.000001) {
      std::cout << "The weights on different tensors must add up to 1!" << std::endl << std::endl ;
      exit(1) ;
    }
  }

  // Initialize the tractography object.
  FilterModel *filter_model = NULL; //Silence warnings.  This will cause segfault if it ever reaches this point.
  Tractography::model_type filter_model_type = Tractography::_1T;

  if (numTensor == 1) {
    if (simpleTensorModel && !freeWater) {
      std::cout << "Using 1-tensor simple model." << std::endl;
      filter_model = new Simple1T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_1T;
    } else if (simpleTensorModel && freeWater) {
      std::cout << "Using 1-tensor simple model with free water estimation." << std::endl;
      filter_model = new Simple1T_FW(Qm, Ql, Qw, Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_1T_FW;
    } else if (!simpleTensorModel && !freeWater) {
      std::cout << "Using 1-tensor full model." << std::endl;
      filter_model = new Full1T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_1T_FULL;
    } else if (!simpleTensorModel && freeWater) {
      std::cout << "Using 1-tensor full model with free water estimation." << std::endl;
      filter_model = new Full1T_FW(Qm, Ql, Qw, Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_1T_FW_FULL;
    }
  } else if (numTensor == 2) {
    if (simpleTensorModel && !freeWater) {
      std::cout << "Using 2-tensor simple model." << std::endl;
      filter_model = new Simple2T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_2T;
    } else if (simpleTensorModel && freeWater) {
      std::cout << "Using 2-tensor simple model with free water estimation." << std::endl;
      filter_model = new Simple2T_FW(Qm, Ql, Qw, Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_2T_FW;
    } else if (!simpleTensorModel && !freeWater) {
      std::cout << "Using 2-tensor full model." << std::endl;
      filter_model = new Full2T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_2T_FULL;
    } else if (!simpleTensorModel && freeWater) {
      std::cout << "Using 2-tensor full model with free water estimation." << std::endl;
      filter_model = new Full2T_FW(Qm, Ql, Qw, Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_2T_FW_FULL;
    }
  } else if (numTensor == 3) {
    if (simpleTensorModel) {
      std::cout << "Using 3-tensor simple model." << std::endl;
      filter_model = new Simple3T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_3T;
    } else {
      std::cout << "Using 3-tensor full model." << std::endl;
      filter_model = new Full3T(Qm, Ql, Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_3T_FULL;
    }
  }



  std::cout << std::endl ;

  Tractography *tract = new Tractography(filter_model, filter_model_type,
                                         tracts, tractsWithSecondTensor,
                                         recordFA, recordNMSE, recordTrace, recordState,
                                         recordCovariance, recordFreeWater, recordTensors,
                                         !noTransformPosition, storeGlyphs, branchesOnly,

                                         minFA, minGA, seedFALimit,
                                         numTensor, seedsPerVoxel,
                                         minBranchingAngle, maxBranchingAngle,
                                         !simpleTensorModel, freeWater,

                                         stepLength, maxHalfFiberLength,
                                         labels,

                                         P0,  SIGMA_SIGNAL, SIGMA_MASK,
                                         MIN_RADIUS, FULL_BRAIN_GA_MIN,

                                         actuallNumThreadsUsed
                                        ) ;

  // if specified on command line, write out binary tract file
  tract->SetWriteBinary(writeBinaryTracts);

  if (tract->LoadFiles(dwiFile, seedsFile, maskFile, normalizedDWIData, outputNormalizedDWIData)) {
    delete tract;
    delete filter_model;
    return 1;
  }

  // Run the tractography.
  tract->Run();

  // Clean up.
  delete tract;
  delete filter_model;

  return 0;
}
