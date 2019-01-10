#include "BRAINSThreadControl.h"

namespace BRAINSUtils
{
  StackPushITKDefaultNumberOfThreads::StackPushITKDefaultNumberOfThreads(const int desiredCount) :
#if ITK_VERSION_MAJOR < 5
    m_originalThreadValue( itk::MultiThreader::GetGlobalDefaultNumberOfThreads() )
#else
    m_originalThreadValue( itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads() )
#endif
  {
    int threadCount(-1);

    if( desiredCount > 0 )  // NOTE: Default is -1, which then uses the ITK default.
      {
      // Respect the users request irregardless of what environmental variables are set to.
      threadCount = desiredCount;
      }
    else
      {                                          // If user set desiredCount <= 0, then use evnironmentanl or internal
                                                 // ITKv4 values.
      threadCount = this->m_originalThreadValue; // This is the old default.
        {                                        // Process the NSLOTS environmental varialble set by the SGE batch
                                                 // processing system
        int         NSLOTSThreadCount(-1);
        std::string numThreads;
        if( itksys::SystemTools::GetEnv("NSLOTS", numThreads) )
          {
          std::istringstream s(numThreads, std::istringstream::in);
          s >> NSLOTSThreadCount;
          }
        if( NSLOTSThreadCount > threadCount )
          {
          threadCount = NSLOTSThreadCount;
          }
        }
      }
    if( threadCount > 0 )
      {
#if ITK_VERSION_MAJOR < 5
      itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threadCount);
#else
      itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(threadCount);
#endif
      }
  }

  StackPushITKDefaultNumberOfThreads::~StackPushITKDefaultNumberOfThreads()
  {
#if ITK_VERSION_MAJOR < 5
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(this->m_originalThreadValue);
#else
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(this->m_originalThreadValue);
#endif
  }
}
