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
#include <vnl/vnl_matrix.h>

/** Calculate fractional anisotropy from eigenvalues */
double l2fa(double l1, double l2, double l3);

/** Calculate Generalized anisotropy from signal */
double s2ga(const vnl_matrix<double>& signal);

/** Calculate curve radius from fiber */
double curve_radius(const std::vector<vec_t>& fiber);

#endif  // UTILITIES_H_
