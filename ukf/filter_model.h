/**
 * \file filter_model.h
 * \brief Contains all of the model reconstructions
 * \todo It would be a minor change to constrain the non free water cases as well. There would be a small additional time cost. Needs to be decided.
*/

#ifndef FILTER_MODEL_H_
#define FILTER_MODEL_H_

#include <cassert>
#include <vector>
#include <iostream>

#include "linalg.h"
#include "ISignalData.h"
#include "unscented_kalman_filter.h"

/**
 * \struct FilterModel
 * \brief Interface for implementation of a signal model
 * \todo IMO rename to SignalModel, FilterModel is confusing.
 * A generic class that defines the transition function and the observation
 * model to be used in the Kalman filter.
*/
class FilterModel
  {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Constructor */
  FilterModel(const int local_state_dim, const ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : _state_dim(local_state_dim),
      _rs(rs), _signal_dim(0), _signal_data(NULL), weights_on_tensors_(weights_on_tensors),
    _constrained(constrained)
  {

    _Q.resize(_state_dim, _state_dim);
    _Q.setConstant(ukfZero); // necessary because otherwise there is memory left overs in the matrix

    if( constrained )
      {
      std::cout << "Using constrained filter\n";
      }
    else
      {
      std::cout << "Using unconstrained filter\n";
      }

  }

  /** Destructor */
  virtual ~FilterModel()
  {
  }

  /** state transition function */
  virtual void F(ukfMatrixType& X) const = 0;

  /** observation, i.e. signal reconstruction */
  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const = 0;

  // Functions that convert the state into a tensor representation that can be
  // used for tractography (meaning the main direction and all the eigenvalues/
  // lambdas).

  /** Extracts principal diffusion direction and eigen values from the state for the 1T cases */
  virtual void State2Tensor1T(const State &, vec3_t &, vec3_t & )
  {
    assert(!"Not implemented");
  }

  /** Extracts principal diffusion direction and eigen values from the state for the 2T cases */
  virtual void State2Tensor2T(const State &, const vec3_t &, vec3_t &,
                            vec3_t &, vec3_t &, vec3_t & )
  {
    assert(!"Not implemented");
  }

  /** Extracts principal diffusion direction and eigen values from the state for the 3T cases */
  virtual void State2Tensor3T(const State &, const vec3_t &, vec3_t &,
                            vec3_t &, vec3_t &, vec3_t &, vec3_t &,
                            vec3_t & )
  {
    assert(!"Not implemented");
  }

  /** Returns the dimension of the state */
  int state_dim() const
  {
    return _state_dim;
  }

  /** Set the dimension of the signal, and resize the corresponding matrices correspondigly */
  void set_signal_dim(const int dim)
  {
    _signal_dim = dim;

    _R.resize(_signal_dim, _signal_dim);
    _R.setConstant(ukfZero); // necessary because otherwise there is memory leftovers in the matrix
    for( int i = 0; i < _signal_dim; ++i )
      {
      _R(i, i) = _rs;
      }

  }

  /** Returns the dimension of the signal */
  int signal_dim() const
  {
    return _signal_dim;
  }

  /** The noise in the state transfer function used by the UKF. */
  const ukfMatrixType & Q() const
  {
    return _Q;
  }

  /** The noise in the signal reconstruction used by the UKF. */
  const ukfMatrixType & R() const
  {
    return _R;
  }

  // The next two functions are only used for the constrained case

  /** The inequality constraint matrix for the constrained UKF */
  const ukfMatrixType & D() const
  {
    return _D;
  }

  /** The inequality constraint right hand side for the constrained UKF */
  const ukfVectorType & d() const
  {
    return _d;
  }

  /** Are we running the contrained version of the filter? */
  bool isConstrained() const
  {
    return _constrained;
  }

  /** Set the pointer to the diffusion signal data */
  void set_signal_data(ISignalData *signal_data)
  {
    _signal_data = signal_data;
  }

protected:

  /** Checks if d is smaller than a small negative threshold. If yes an error is returned. Otherwise d is rounded to ukfZero
    */
  ukfPrecisionType CheckZero(const ukfPrecisionType & d) const;

  /** The dimension of the state */
  const int _state_dim;

  /** The constant signal noise for each signal component */
  ukfPrecisionType _rs;

  /** The dimension of the signal */
  int _signal_dim;

  /** Pointer to the diffusion signal data */
  ISignalData *_signal_data;

  /** Process noise */
  ukfMatrixType _Q;

  /** Signal reconstruction noise */
  ukfMatrixType _R;

  /** Inequality constraint matrix, only used for constrained UKF */
  ukfMatrixType _D;

  /** Inequality right hand side, only used for constrained UKF */
  ukfVectorType _d;

  /** The weights of each tensor. (Most commonly equal all are equal) */
  const ukfVectorType weights_on_tensors_;

  /** Are we using the constrained filter */
  bool _constrained;
  };

inline mat33_t SetIdentityScaled(double diff_fw)
{
   mat33_t tmp;
   tmp.setIdentity();
   tmp*= diff_fw;
   return tmp;
}

/* Common Function used for NODDI */
extern void createProtocol(const ukfVectorType& _b_values,
                    ukfVectorType& _gradientStrength, ukfVectorType& _pulseSeparation);

#endif  // FILTER_MODEL_H_
