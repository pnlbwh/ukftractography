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
double l2fa(double l1, double l2, double l3);

/** Calculate Generalized anisotropy from signal */
double s2ga(const ukfMatrixType& signal);

/** Calculate curve radius from fiber */
double curve_radius(const stdVec_t& fiber);

#endif  // UTILITIES_H_
