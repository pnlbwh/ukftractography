/**
 * \file utilities.h
 * \brief Calculation of frequently used (fa,ga,curve_radius) defined in this file
*/
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cassert>
#include <cmath>
#include <vector>
#include "linalg.h"

/** Calculate fractional anisotropy from eigenvalues */
ukfPrecisionType l2fa(ukfPrecisionType l1, ukfPrecisionType l2, ukfPrecisionType l3);

/** Calculate Generalized anisotropy from signal */
ukfPrecisionType s2ga(const ukfMatrixType& signal);

/** Calculate curve radius from fiber */
ukfPrecisionType curve_radius(const stdVec_t& fiber, int size);

#endif  // UTILITIES_H_
