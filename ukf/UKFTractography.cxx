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

#include "cli.h"
#include "tractography.h"
#include "BRAINSThreadControl.h"

#include "vtkNew.h"
#include "vtkPolyData.h"

extern "C" {

UKFTRACTOGRAPHYLIB_EXPORT int ModuleEntryPoint(int argc, char **argv)
{
  UKFSettings ukf_settings;

  if (int stat = ukf_parse_cli(argc, argv, ukf_settings) != EXIT_SUCCESS)
    return stat;

  // NOTE:  When used as share libary one must be careful not to permanently reset number of threads
  //        for entire program (i.e. when used as a slicer modules.
  //        This also addresses the issue when the program is run as part of a batch processing
  //        system so that the number of cores allocated by scheduler is respected rather than
  //        blindly using all the cores that are found.
  //        This implementation is taken from extensive testing of the BRAINSTools
  // (this object will be deleted by RAII and return to original thread count)
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(ukf_settings.num_threads);
  const int actualNumThreadsUsed = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
  ukf_settings.num_threads = actualNumThreadsUsed;
  {
    std::cout << "Found " << actualNumThreadsUsed << " cores on your system." << std::endl;
    std::cout << "Running tractography with " << actualNumThreadsUsed << " thread(s)." << std::endl;
  }

  // TODO these have always been hard-coded here. But why?
  bool normalizedDWIData = false;
  bool outputNormalizedDWIData = false;

  // initializing super object
  Tractography *tract = new Tractography(ukf_settings);

  // if specified on command line, write out binary tract file
  tract->SetWriteBinary(!ukf_settings.writeAsciiTracts);
  tract->SetWriteCompressed(!ukf_settings.writeUncompressedTracts);

  int writeStatus = 0;
  try
    {
    if (tract->LoadFiles(ukf_settings.dwiFile,
                         ukf_settings.seedsFile,
                         ukf_settings.maskFile,
                         normalizedDWIData, outputNormalizedDWIData) == EXIT_FAILURE)
      {
      itkGenericExceptionMacro(<< "::LoadFiles failed with unknown error.");
      }

    // set filter model TODO: refactor?
    tract->SetFilterModelType((Tractography::model_type)ukf_settings.filter_model_type);

    // Run the tractography.
    writeStatus = tract->Run();
    std::cout << "H count = " << countH << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "UKFTractography: ITK ExceptionObject caught!" << std::endl;
    std::cerr << err << std::endl;

    writeStatus = EXIT_FAILURE;
    }
  catch( std::exception& exc )
    {
    std::cerr << "UKFTractography: std::exception caught:" << std::endl;
    std::cerr << exc.what() << std::endl;

    writeStatus = EXIT_FAILURE;
    }
  catch(...)
    {
    std::cerr << "UKFTractography: Unknown exception caught!" << std::endl;

    writeStatus = EXIT_FAILURE;
    }

  // Clean up.
  delete tract;

  return writeStatus;
}

}; // extern "C"
