/**
 * \file utilities.h
 * \brief Calculation of frequently used (fa,ga,curve_radius) defined in this file
 * Calculations for NODDI model are implemented
*/
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <cassert>
#include <cmath>
#include <vector>
#include "linalg.h"

#define NaN std::numeric_limits<double>::quiet_NaN()

/** Calculate fractional anisotropy from eigenvalues */
ukfPrecisionType l2fa(ukfPrecisionType l1, ukfPrecisionType l2, ukfPrecisionType l3);

/** Calculate Generalized anisotropy from signal */
ukfPrecisionType s2ga(const ukfMatrixType& signal);

/** Calculate Generalized anisotropy from signal */
ukfPrecisionType s2adc(const ukfMatrixType& signal);

/** Calculate curve radius from fiber */
ukfPrecisionType curve_radius(const stdVec_t& fiber);

// special case for real x
double dawsonf(double kappa);

/* returns the equivalent parallel and perpendicular diffusion coefficients
   for hindered compartment with impermeable cylinder's oriented with a
   Watson's distribution with a cocentration parameter of kappa */
void WatsonHinderedDiffusionCoeff(ukfPrecisionType dPar, ukfPrecisionType dPerp, ukfPrecisionType kappa, ukfMatrixType& dw);

void ExtraCelluarModel(ukfPrecisionType dPar, ukfPrecisionType Vic, ukfPrecisionType kappa,
                       ukfVectorType& gradientStrength, ukfVectorType& pulseSeparation,
                       const stdVec_t& u, vec3_t& fiberdir, ukfVectorType& Eec);

void IntraCelluarModel(ukfPrecisionType dPar, ukfPrecisionType kappa, ukfVectorType& gradientStrength,
                       ukfVectorType& pulseSeparation, const stdVec_t& u, vec3_t& fiberdir, ukfVectorType& Eic);

void IsoModel(ukfPrecisionType dIso, ukfVectorType& gradientStrength, ukfVectorType& pulseSeparation,
              ukfVectorType& Eiso);

#endif  // UTILITIES_H_
