/**
 * \file UKFTractography.cxx
 * \brief Main file of the project
 *
 * In this file the input arguments are parsed and passed on to the tractography
 * Also the model choice happens here
*/

unsigned int countH=0;
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "BRAINSThreadControl.h"
#include "filter_Full1T.h"
#include "filter_Full1T_FW.h"
#include "filter_Full2T.h"
#include "filter_Full2T_FW.h"
#include "filter_Full3T.h"
#include "filter_NODDI1F.h"
#include "filter_NODDI2F.h"
#include "filter_Simple1T.h"
#include "filter_Simple1T_FW.h"
#include "filter_Simple2T.h"
#include "filter_Simple2T_FW.h"
#include "filter_Simple3T.h"
#include "filter_Simple2BiExp_FW.h"
#include "tractography.h"
#include "UKFTractographyCLP.h"


static const bool verbose = true;

void setAndTell(ukfPrecisionType & x, const ukfPrecisionType y, const std::string & name)
{
  if (verbose) {
    x = y;
    std::cout << "- " << name << ": " << y << std::endl;
  }
}

void tell(const ukfPrecisionType & x, const std::string &name)
{
  if (verbose) {
    std::cout << "* " << name << ": " << x << std::endl;
  }
}

int main(int argc, char **argv)
{

  PARSE_ARGS ;
  ukfPrecisionType l_minFA = minFA;
  ukfPrecisionType l_minGA = minGA;
  ukfPrecisionType l_stepLength = stepLength;
  ukfPrecisionType l_recordLength = recordLength;
  ukfPrecisionType l_maxHalfFiberLength = maxHalfFiberLength;
  ukfPrecisionType l_seedFALimit = seedFALimit;
  ukfPrecisionType l_Qm = Qm;
  ukfPrecisionType l_Ql = Ql;
  ukfPrecisionType l_Qt = Qt;
  ukfPrecisionType l_Qw = Qw;
  ukfPrecisionType l_Qkappa = Qkappa;
  ukfPrecisionType l_Qvic = Qvic;
  ukfPrecisionType l_Rs = Rs;
  ukfPrecisionType l_maxBranchingAngle = maxBranchingAngle;
  ukfPrecisionType l_minBranchingAngle = minBranchingAngle;

  // If sigmaSignal is not set minimum of voxel size is used for interpolation
  ukfPrecisionType SIGMA_SIGNAL = sigmaSignal;

  bool simpleTensorModel = !fullTensorModel;


  std::cout << std::endl;

  // CONSTANTS
  const ukfPrecisionType SIGMA_MASK 	      = 0.5;
  const ukfPrecisionType P0 			          = 0.01;
  const ukfPrecisionType MIN_RADIUS 		    = 0.87;
  const ukfPrecisionType FULL_BRAIN_GA_MIN	= 0.18;
  const ukfPrecisionType D_ISO			        = 0.003; // Diffusion coefficient of free water

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

  if (l_maxHalfFiberLength <= 0) {
    std::cout << "Invalid maximum half fiber length!" << std::endl ;
    return 1 ;
  }

  if (noddi) {
    if (recordFA || recordTrace || recordFreeWater || recordTensors) {
        std::cout << "recordFA, recordTrace, recordFreeWater, recordTensors parameters can only be used with tensor models\n";
        return 1 ;
    }
  }
  if (diffusionPropagator) { 
    if (noddi) {
      std::cout<<"Noddi and Diffusion Propagator parameters are mutually exclusive. Use either one of the two"<<std::endl;
      return 1;
    } 
    
    if (recordFA || recordTrace) {
      std::cout << "recordFA and recordTrace cannot be used with the diffusion propagator model\n";
      return 1 ;
    }
     
    if (!freeWater) {
      std::cout<<"Since the Diffusion Propagator model is used, the free water parameter will be estimated"<<std::endl;
      freeWater = true;
    }
    
    if (!simpleTensorModel) {
      std::cout<<"Since the Diffusion Propagator model is used, the simple tensor model will be used"<<std::endl;
      simpleTensorModel = true;
    }
    
    if (numTensor != 2) {
      std::cout<<"Since the Diffusion Propagator model is used, the number of tensors is set to two (2)"<<std::endl;
      numTensor = 2;
    }    
  
  }
  
  if (recordRTOP && !diffusionPropagator) {
    std::cout<<"recordRTOP cannot be used with any other models than the diffusionPropagator"<<std::endl;
    return 1;
  }
  
  if (l_Qt != 0.0 && !diffusionPropagator) {
    std::cout<<"Qt parameter cannot be set with any other models than the diffusionPropagator model"<<std::endl;
    return 1; 
  }
  
  if (minRTOP != 0.0 && !diffusionPropagator) {
    std::cout<<"minRTOP parameter cannot be set with any other models than the diffusionPropagator model"<<std::endl;
    return 1; 
  }
  
  if (maxNMSE != 0.0 && !diffusionPropagator) {
    std::cout<<"maxNMSE parameter cannot be set with any other models than the diffusionPropagator model"<<std::endl;
    return 1;   
  }
  
  if (maxUKFIterations != -1.0 && !diffusionPropagator) {
    std::cout<<"maxUKFIterations parameter cannot be set with any other models than the diffusionPropagator model"<<std::endl;
    return 1;   
  }
  
  
//   if (l_stepLength <= 0){
//     std::cout << "Invalid step length!" << std::endl ;
//     return 1 ;
//   }

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
  std::string strModel = simpleTensorModel ? "simple model" : "full model";
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

  if (l_minFA == 0.15) {
    setAndTell(l_minFA, l_minFA, "minFA");
  } else {
    tell(l_minFA, "minFA");
  }

  if (l_seedFALimit == 0.0) {
    setAndTell(l_seedFALimit, FULL_BRAIN_GA_MIN, "seedFALimit");  // Used to default to 2 times the FA threshold.
  } else {
    tell(l_seedFALimit, "seedFALimit");
  }

  if(!noddi){
    if(recordVic || recordKappa || recordViso){
      std::cout << "Can use recordVic or recordKappa or recordViso parameters only with noddi model";
      exit(1);
    }
  }

  if (l_Qm == 0.0) {
    if (noddi){
      if (numTensor == 1)
        setAndTell(l_Qm, 0.0025, "Qm");
      else
        setAndTell(l_Qm, 0.001, "Qm");
    } 
    else if (diffusionPropagator) {
        setAndTell(l_Qm, 0.0001, "Qm");
    }
    else if (numTensor == 1) {
        setAndTell(l_Qm, 0.005, "Qm");//l_Qm = 0.0015;
    } else {
      if (!simpleTensorModel) {
        setAndTell(l_Qm, 0.002, "Qm");//l_Qm = 0.002;
      } else {
        setAndTell(l_Qm, 0.001, "Qm");//l_Qm = 0.001; was 0.003, changed to 0.001 for new Interp3Signal
      }
    }
  } else {
    tell(l_Qm, "Qm");
  }


  if (noddi){
    if ( l_Qkappa == 0.0)
      setAndTell (l_Qkappa, 0.01, "Qkappa");
    else
      tell(l_Qkappa, "Qkappa");
  }
  else {
    if (l_Ql == 0.0) {
      if (diffusionPropagator) {
        setAndTell(l_Ql, 150, "Ql");
      }
      else if (numTensor == 1) {
        setAndTell(l_Ql, 300.0, "Ql");//l_Ql = 25.0;
      } else if (numTensor == 2) {
        setAndTell(l_Ql, 50.0, "Ql");//was l_Ql = 100.0; for old Interp3Signal
      } else if (numTensor == 3) {
        setAndTell(l_Ql, 100.0, "Ql");//l_Ql = 150.0;
      }
    } else {
        tell(l_Ql, "Ql");
    }
  }
 // In the diffusion propagator model, Ql is used for the fast diffusion eigenvalues, and Qt for the slow diffusion.
  if (diffusionPropagator) {
    if (l_Qt == 0.0) {
      setAndTell(l_Qt, 50, "Qt");
    }
    else {
      tell(l_Qt, "Qt");
    }
  }
  
  if (diffusionPropagator) {
    if (minRTOP == 0.0) {
      setAndTell(minRTOP, 60.0, "minRTOP");
    }
    else {
      tell(minRTOP, "minRTOP");
    }
  }
  
  if (diffusionPropagator) {
    if (maxNMSE == 0.0) {
      setAndTell(maxNMSE, 0.15, "maxNMSE");
    }
    else {
      tell(maxNMSE, "maxNMSE");
    }
  }
  
  if (diffusionPropagator) {
    if (maxUKFIterations == -1.0) {
      setAndTell(maxUKFIterations, 5, "maxUKFIterations");
    }
    else {
    if (maxUKFIterations < 0.0) {
        std::cout<<"Error: maxUKFIterations cannot be negative. Exiting"<<std::endl;
        exit(1);
    }
      tell(maxUKFIterations, "maxUKFIterations");
    }
  }  


  if (l_Rs == 0.0) {
    if (diffusionPropagator) {
      setAndTell(l_Rs, 0.015, "Rs");
    }
    else if (numTensor == 1) {
      setAndTell(l_Rs, 0.01, "Rs");//l_Rs = 0.02;
    } else {
      if (!simpleTensorModel) {
        setAndTell(l_Rs, 0.01, "Rs");// = 0.01;
      } else {
        setAndTell(l_Rs, 0.02, "Rs");//was l_Rs = 0.015;for old Interp3Signal
      }
    }
  } else {
    tell(l_Rs, "Rs");
  }

  if (l_stepLength == 0.3 && diffusionPropagator)  {
    setAndTell(l_stepLength, 0.5, "stepLength");
  }
  else if (l_stepLength == 0.0) {
    if (diffusionPropagator) {
      setAndTell(l_stepLength, 0.5, "stepLength");
    }
    else if (numTensor == 1) {
      setAndTell(l_stepLength, 0.3, "stepLength");
    } else if (numTensor == 2) {
      setAndTell(l_stepLength, 0.3, "stepLength"); //was 0.2 for old Interp3Signal
    } else { // 3T
      setAndTell(l_stepLength, 0.15, "stepLength");
    }
  } else {
    tell(l_stepLength, "stepLength");
  }

  if (l_recordLength == 0.0) {
    if (numTensor == 1) {
      setAndTell(l_recordLength, 0.9, "recordLength");
    } else if (numTensor == 2) {
      setAndTell(l_recordLength, 0.9, "recordLength"); //was 0.2 for old Interp3Signal
    } else { // 3T
      setAndTell(l_recordLength, 0.45, "recordLength");
    }
  } else {
    tell(l_recordLength, "recordLength");
  }

  if (noddi){
    if (l_Qvic == 0.0)
      if (numTensor == 1)
        setAndTell(l_Qvic, 0.0005, "Qvic = Qviso");
      else
        setAndTell(l_Qvic, 0.004, "Qvic = Qviso");
    else
      tell(l_Qvic, "Qvic = Qviso");
  }
  else
    if (freeWater) {
      if (l_Qw == 0.0) {
        if (diffusionPropagator) {
          setAndTell(l_Qw, 0.002, "Qw");
        }
        else if (numTensor == 1) {
          setAndTell(l_Qw, 0.0025, "Qw"); // estimated in a paramsearch // 0.0025
        } else if (numTensor == 2) {
          setAndTell(l_Qw, 0.0015, "Qw"); // 0.0015
        }
      } else {
          tell(l_Qw, "Qw");
      }
    }

  tell(l_minGA, "minGA");

  if (seedsPerVoxel == 1) {
    std::cout << "- seedsPerVoxel: " << seedsPerVoxel << std::endl;
  } else {
    std::cout << "* seedsPerVoxel: " << seedsPerVoxel << std::endl;
  }
  bool noTransformPosition = false;
  bool branchesOnly = false;
  //if (normalizedDWIData) {
    //outputNormalizedDWIData = false ;
  //}
  bool normalizedDWIData = false;
  bool outputNormalizedDWIData = false;

  ukfVectorType weightsOnTensors(numTensor);
  for (int i = 0; i < numTensor; i++)
    {
    weightsOnTensors[i]=(1.0 / numTensor) ;
    }

  ukfPrecisionType weight_accumu = 0 ;
  for (int i = 0; i < numTensor; i++)
    {
    weight_accumu += weightsOnTensors[i] ;
    }
  if (std::abs(weight_accumu - 1.0) > 0.000001)
    {
    std::cout << "The weights on different tensors must add up to 1!" << std::endl << std::endl ;
    exit(1) ;
    }
  else
    {
    weightsOnTensors.norm(); // Normalize for all to add up to 1.
    }


  // Initialize the tractography object.
  FilterModel *filter_model = NULL; //Silence warnings.  This will cause segfault if it ever reaches this point.
  Tractography::model_type filter_model_type = Tractography::_1T;
  
  if (diffusionPropagator) {
    std::cout << "Using Diffusion Propagator model (2 simple tensors, free water, bi-exponential)" << std::endl;
    filter_model = new Simple2T_BiExp_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO, l_Qt);
    filter_model_type = Tractography::_2T_BiExp_FW;
    
  } else if (noddi){
    if (numTensor == 1){
      std::cout << "Using NODDI 1-Fiber model." << std::endl;
      filter_model = new NODDI1F(l_Qm, l_Qkappa, l_Qvic, l_Rs, weightsOnTensors, noddi);
      filter_model_type = Tractography::_1T_FW; // same vtk writer can be used
    } else if (numTensor == 2){
      std::cout << "Using NODDI 2-Fiber model." << std::endl;
      filter_model = new NODDI2F(l_Qm, l_Qkappa, l_Qvic, l_Rs, weightsOnTensors, noddi);
      filter_model_type = Tractography::_2T_FW; // same vtk writer can be used
    }
  }else if (numTensor == 1) {
    if (simpleTensorModel && !freeWater) {
      std::cout << "Using 1-tensor simple model." << std::endl;
      filter_model = new Simple1T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_1T;
    } else if (simpleTensorModel && freeWater) {
      std::cout << "Using 1-tensor simple model with free water estimation." << std::endl;
      filter_model = new Simple1T_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_1T_FW;
    } else if (!simpleTensorModel && !freeWater) {
      std::cout << "Using 1-tensor full model." << std::endl;
      filter_model = new Full1T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_1T_FULL;
    } else if (!simpleTensorModel && freeWater) {
      std::cout << "Using 1-tensor full model with free water estimation." << std::endl;
      filter_model = new Full1T_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_1T_FW_FULL;
    }
  } else if (numTensor == 2) {
    if (simpleTensorModel && !freeWater) {
      std::cout << "Using 2-tensor simple model." << std::endl;
      filter_model = new Simple2T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_2T;
    } else if (simpleTensorModel && freeWater) {
      std::cout << "Using 2-tensor simple model with free water estimation." << std::endl;
      filter_model = new Simple2T_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_2T_FW;
    } else if (!simpleTensorModel && !freeWater) {
      std::cout << "Using 2-tensor full model." << std::endl;
      filter_model = new Full2T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_2T_FULL;
    } else if (!simpleTensorModel && freeWater) {
      std::cout << "Using 2-tensor full model with free water estimation." << std::endl;
      filter_model = new Full2T_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO);
      filter_model_type = Tractography::_2T_FW_FULL;
    }
  } else if (numTensor == 3) {
    if (simpleTensorModel) {
      std::cout << "Using 3-tensor simple model." << std::endl;
      filter_model = new Simple3T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_3T;
    } else {
      std::cout << "Using 3-tensor full model." << std::endl;
      filter_model = new Full3T(l_Qm, l_Ql, l_Rs, weightsOnTensors, freeWater);
      filter_model_type = Tractography::_3T_FULL;
    }
  }



  std::cout << std::endl ;

  Tractography *tract = new Tractography(filter_model, filter_model_type,
                                         tracts, tractsWithSecondTensor,
                                         recordFA, recordNMSE, recordTrace, recordState,
                                         recordCovariance, recordFreeWater, recordTensors,
                                         recordVic, recordKappa, recordViso,
                                         !noTransformPosition, storeGlyphs, branchesOnly,

                                         l_minFA, l_minGA, l_seedFALimit,
                                         numTensor, seedsPerVoxel,
                                         l_minBranchingAngle, l_maxBranchingAngle,
                                         !simpleTensorModel, freeWater, noddi,
                                         diffusionPropagator, minRTOP, recordRTOP,
                                         maxNMSE, maxUKFIterations,

                                         l_stepLength, l_recordLength, l_maxHalfFiberLength,
                                         labels,

                                         P0,  SIGMA_SIGNAL, SIGMA_MASK,
                                         MIN_RADIUS, FULL_BRAIN_GA_MIN,

                                         actuallNumThreadsUsed
                                        ) ;

  // if specified on command line, write out binary tract file
  tract->SetWriteBinary(!writeAsciiTracts);
  tract->SetWriteCompressed(!writeUncompressedTracts);

  if (tract->LoadFiles(dwiFile, seedsFile, maskFile, normalizedDWIData, outputNormalizedDWIData)) {
    delete tract;
    delete filter_model;
    return 1;
  }

  // Run the tractography.
  const int writeStatus = tract->Run();

  // Clean up.
  delete tract;
  delete filter_model;

  std::cout << "H count = " << countH << std::endl;
  return writeStatus;
}
