/**
 * \file unscented_kalman_filter.h
 * \brief The implementation of the unscented Kalman Filter
 * The C++ implementation of the unscented Kalman Filter (UKF) with the option to run
 * with constraints
*/

#ifndef UNSCENTED_KALMAN_FILTER_H_
#define UNSCENTED_KALMAN_FILTER_H_

#include <vector>
#include "ukf_types.h"

class FilterModel;

/**
 * \class UnscentedKalmanFilter
 * \brief The C++ implementation of the unscented Kalman Filter
*/
class UnscentedKalmanFilter
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * \brief Constructor
   * \param filter_model The UKF is initialized with a filter model that defines the observation
   * and transition functions. Also required are the noise parameters.
  */
  UnscentedKalmanFilter(FilterModel *filter_model);

  /**
   * \brief Does the one filter step.
   * \param[in]  x Current state vector
   * \param[in]  p Current convariance matrix of the stateal and it
   * \param[in]  z The signal interpolated at the current position
   * \param[out] x_new Updated state
   * \param[out] p_new Updated covariance
   * \param[out] The normalized mean squared reconstruction error
  */
  void Filter(const State& x, const ukfMatrixType& p, const ukfVectorType& z, // This is the signal
              State& x_new, ukfMatrixType& p_new, ukfPrecisionType& dNormMSE);

private:
  /** Spreads the points around the current state using the covariance. */
  void SigmaPoints(const State& x, const ukfMatrixType& p, ukfMatrixType& x_spread);

  /**
   * \brief Contrains the state matrix
   * \param X The state matrix which is the result of spreading out the original state through SigmaPoints.
   *          Will be constrained in this function.
   * \param W The covariance necesseray for contraining. See the malcolm MICCAI paper.
  */
  void Constrain(ukfMatrixType& X, const ukfMatrixType& W);

  /**
   * \brief Contrains the state vector
   * \param X The state vector which will be constrained.
   * \param W The covariance necesseray for contraining. See the malcolm MICCAI paper.
  */
  void Constrain(ukfVectorType& x, const ukfMatrixType& W);

  /** A helper function to check if the constrain operation is necessary */
  bool violatesContraints(ukfVectorType& x);

  /** Pointer to the filter model */
  const FilterModel * const m_FilterModel;

  /** state vector dimension */
  // HACK REMOVE int _state_dim;

  /** for the distribution of the sigma ponts (:= sqrt(dim + m_SigmaPointSpread) ) */
  ukfPrecisionType m_Scale;

  /** The weights for spreading the sigma points */
  ukfVectorType m_Weights;

  /** Matrix of weights for spreading of sigma points consisting of the repeted entries of m_Weights */
  ukfMatrixType m_WeightsRepeated;

  /** A fixed parameters used for spreading of the sigma points */
  ukfPrecisionType m_SigmaPointSpread;

  //DUMMY VARIABLES TAHT ARE ALWAYS ZERO and not used.
  ukfMatrixType m_DummyZeroCE;
  ukfVectorType m_DummyZeroce0;
};

#endif  // UNSCENTED_KALMAN_FILTER_H_
