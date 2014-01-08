#ifndef __computdispersion_h
#define __computdispersion_h
#include "fiberbundle.h"

int
computedispersion(fiberbundle &bundle,
                  double scale,
                  unsigned int numberOfSamplingDirections,
                  const std::string &outputFilename,
                  unsigned int tractSubSampling = 1,
                  unsigned int fiberPointSubSampling = 1);

#endif // __computdispersion_h

