/**
 * \file fiber.h
 * \brief Description of a fiber
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#ifndef UKFFIBER_H_
#define UKFFIBER_H_

#include <vector>
#include <cassert>
#include "unscented_kalman_filter.h"
#include "linalg.h"

#include <vnl/vnl_matrix.h>

/**
 * \struct UKFFiber
 * \brief Points of a fiber, and scalars corresponding to the points
 *
 * The points that make a fiber are defined as a vector of 3D points. In addition
 * there is a vector of scalar values for each scalar value that can be recorded
 * of the same length
*/
struct UKFFiber
  {

  /** vector of 3D points defining the fiber path */
  std::vector<vec_t> position;
  /** FA of tensor 1 */
  std::vector<double> fa;
  /** FA of tensor 2 */
  std::vector<double> fa2;
  /** Array 2 norm of the covariance matrix */
  std::vector<double> norm;
  /** State of the current model at the current position*/
  std::vector<State> state;
  /** dim(state) x dim(state) matrix */
  std::vector<ukfMatrixType > covariance;
  /** Percentage of free water i.e. 1-w */
  std::vector<double> free_water;
  /** Normalized mean squared error of the signal reconstruction to the signal */
  std::vector<double> normMSE;
  /** Trace of tensor 1 */
  std::vector<double> trace;
  /** Trace of tensor 2 */
  std::vector<double> trace2;
  };

/**
 * \struct BranchingSeedAffiliation
 * \brief Which fibers belong together
 *
 * The structure to document on which primary fiber and at which position is the seed of the branch
*/
struct BranchingSeedAffiliation
  {
  size_t fiber_index_;
  int position_on_fiber_;
  };

/**
 * \brief Joins two fibers originating from the same seed point
 *
 * A pair of two primary fibers are started from each seed point in two opposite directions. This functions joins them up pairly to
 * form complete primary fibers, and eliminates fibers that are too short. Besides, each branch is back traced to form a whole fiber
*/
void PostProcessFibers( const std::vector<UKFFiber>& raw_primary, const std::vector<UKFFiber>& raw_branch,
                        const std::vector<BranchingSeedAffiliation>& branching_seed_affiliation,
                        const bool branches_only, std::vector<UKFFiber>& fibers);

/** The minimum number of points on a fiber. UKFFiber with fewer points are rejected */
const int MINIMUM_NUM_POINTS_ON_FIBER = 10;

/** Used to get rid of branches, that start near the end of primary fibers. See fiber.cc:70. */
const int FIBER_TAIL_THRESHOLD = 5;

#endif
