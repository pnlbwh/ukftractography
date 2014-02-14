/** Utilities to compute Spherical-Harmonics (SH)-related functions. Take a look to the
 corresponding implementation (.cxx) file for details on the meaning and usage of each
 of the functions in this file.

 IMPORTANT NOTE: Part of the implementations of the functions in this headers file is
 in the file sh2hot.cxx. Hence, it is mandatory to include sh2hot.cxx in any project
 including "sphericalHarmonics.h"
 */

#ifndef __SphHarm_h
#define __SphHarm_h
#include "CompressedSensing.h"

namespace shmaths
{

void computeSHMatrix( const unsigned int N,
                      const double* theta,
                      const double* phi,
                      const unsigned int L,
                      MatrixType& sh );

void computeSHMatrixSymmetric( const unsigned int N,
                               const double* theta,
                               const double* phi,
                               const unsigned int L,
                               MatrixType& sh );

void computeSHEigMatrix( const unsigned int L, MatrixType& eig );

void computeSHEigMatrixSymmetric( const unsigned int L, MatrixType& eig );

void computeSHFRTMatrix( const unsigned int L, MatrixType& frt);

void computeSHFRTMatrixSymmetric( const unsigned int L, MatrixType& frt );

void computeSphericalCoordsFromCartesian( const double x,
                                          const double y,
                                          const double z,
                                          double& theta,
                                          double& phi );

void computeSphericalCoordsFromCartesian( const double* x,
                                          const double* y,
                                          const double* z,
                                          double* theta,
                                          double* phi,
                                          const unsigned int N);
} // End namespace shmaths

#endif //__SphHarm_h
