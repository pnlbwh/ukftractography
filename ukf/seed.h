/**
 * \file seed.h
 * \brief Contains data structure for storing seed point information.
*/

#ifndef SEED_H_
#define SEED_H_

#include "unscented_kalman_filter.h"
#include "linalg.h"

#include <vnl/vnl_matrix.h>

/**
 * \struct SeedPointInfo
 * \brief Describes a seed point
 *
 * Stores all information for a seed point to start tractography. The
 * start_dir is only used for the simple model and the angles only for the
 * complex/full representation.
*/
struct SeedPointInfo {
  /** The state of the state-space represenation of the model. */
  State state;
  /** The covariance matrix of the state */
  vnl_matrix<double> covariance;
  /** The location of the seed */
  vec_t point;
  /** The starting direction for the simple model */
  vec_t start_dir;
  /** Fractional Anisotropy of the first tensor */
  double fa;
  /** Fractional Anisotropy of the second tensor */
  double fa2;
  /** Trace of the first tensor */
  double trace;
  /** Trace of the second tensor */
  double trace2;
} ;

/** Writes debug information about seeds to stdout. */
void PrintSeedInfo(const std::vector<SeedPointInfo>& seed_infos) ;

#endif
