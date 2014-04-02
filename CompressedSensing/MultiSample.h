#ifndef __MultiSample_h
#define  __MultiSample_h
#include "CompressedSensing.h"
#include <vector>
extern void MultiSample(unsigned J,
                        std::vector<MatrixType> &v,
                        unsigned long &M,
                        int m0 = 4);
#endif // __MultiSample_h

