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

/**
 * \struct Full1T_FW
 * \brief Full 1-Tensor model with free water
 *
 * Model describing 1-tensor tractography with the full tensor representation (3 angles, 3 eigenvalues)
 * and free water estimation.
*/
class Full1T_FW : public FilterModel
  {
public:
  Full1T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained,
            const ukfPrecisionType diff_fw)
    : FilterModel(7, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw))
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
    _Q(6, 6) = qw;

    _D.resize(7, 5);
    _D.setConstant(ukfZero);

    _d.resize(5);

    // Setting the constraints according to D'*x >= -d
    _D(6, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(6, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2) = 1;
    _d(2) = 0;  // l1 >= 0
    _D(4, 3) = 1;
    _d(3) = 0;  // l2 >= 0
    _D(5, 4) = 1;
    _d(4) = 0;  // l3 >= 0

  }

  virtual ~Full1T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;

  };

/**
 * \struct Full1T
 * \brief Full 1-Tensor model without free water
 *
 * Model describing 1-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
class Full1T : public FilterModel
  {
public:
  Full1T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
  }

  virtual ~Full1T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;
  };

/**
 * \struct Full2T
 * \brief Full 2-Tensor model without free water
 *
 * Model describing 2-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
class Full2T : public FilterModel
  {
public:
  Full2T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(12, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
  }

  virtual ~Full2T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;
  };

/**
 * \struct Full2T_FW
 * \brief Full 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the full tensor representation (3 angles, 3 eigenvalues)
 * and free water estimation.
*/
class Full2T_FW : public FilterModel
  {
public:
  Full2T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained,
            const ukfPrecisionType diff_fw)
    : FilterModel(13, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw))
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
    _Q(12, 12) = qw; // noise for weights

    _D.resize(13, 8);
    _D.setConstant(ukfZero);

    _d.resize(8);

    // Setting the constraints according to D'*x >= -d
    _D(12, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(12, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2)  = 1;
    _d(2) = 0;  // l11 >= 0
    _D(4, 3)  = 1;
    _d(3) = 0;  // l21 >= 0
    _D(5, 4)  = 1;
    _d(4) = 0;  // l31 >= 0
    _D(9, 5)  = 1;
    _d(5) = 0;  // l12 >= 0
    _D(10, 6) = 1;
    _d(6) = 0;  // l22 >= 0
    _D(11, 7) = 1;
    _d(7) = 0;  // l32 >= 0
  }

  virtual ~Full2T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  const mat33_t m_D_iso;
  };

/**
 * \struct Full3T
 * \brief Full 3-Tensor model with free water
 *
 * Model describing 3-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
class Full3T : public FilterModel
  {
public:
  Full3T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(18, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(12, 12) = _Q(13, 13) = _Q(14, 14) = qs;

    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
    _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
    _Q(15, 15) = _Q(16, 16) = _Q(17, 17) = ql;
  }

  virtual ~Full3T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor3T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2, vec3_t& m3,
                            vec3_t& l3);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  };

/**
 * \struct Simple1T_FW
 * \brief Simple 1-Tensor model with free water
 *
 * Model describing 1-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
class Simple1T_FW : public FilterModel
  {
public:
  Simple1T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
              bool constrained, const ukfPrecisionType diff_fw)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw))
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = ql;
    _Q(5, 5) = qw; // noise for weights

    _D.resize(6, 4);
    _D.setConstant(ukfZero);

    _d.resize(4);

    // Setting the constraints according to D'*x >= -d
    _D(5, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(5, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2) = 1;
    _d(2) = 0;  // l1 >= 0
    _D(4, 3) = 1;
    _d(3) = 0;  // l2 >= 0

  }

  virtual ~Simple1T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;

  };

/**
 * \struct Simple1T
 * \brief Simple 1-Tensor model
 *
 * Model describing 1-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
*/
class Simple1T : public FilterModel
  {
public:
  Simple1T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(5, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = ql;
  }

  virtual ~Simple1T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  const ukfPrecisionType _lambda_min;
  };

/**
 * \struct Simple2T
 * \brief Simple 2-Tensor model
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal).
*/
class Simple2T : public FilterModel
  {
public:
  Simple2T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(10, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
  }

  virtual ~Simple2T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  };

/**
 * \struct Simple2T_FW
 * \brief Simple 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
class Simple2T_FW : public FilterModel
  {
public:
  Simple2T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
              bool constrained, const ukfPrecisionType diff_fw)
    : FilterModel(11, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw))
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
    _Q(10, 10) = qw; // noise for weights

    _D.resize(11, 6);
    _D.setConstant(ukfZero);

    _d.resize(6);

    // Setting the constraints according to D'*x >= -d
    _D(10, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(10, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2)  = 1;
    _d(2) = 0;  // l11 >= 0
    _D(4, 3)  = 1;
    _d(3) = 0;  // l21 >= 0
    _D(8, 4)  = 1;
    _d(4) = 0;  // l12 >= 0
    _D(9, 5)  = 1;
    _d(5) = 0;  // l22 >= 0
  }

  virtual ~Simple2T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;
  };

/**
 * \struct Simple3T
 * \brief Simple 3-Tensor model
 *
 * Model describing 3-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
*/
class Simple3T : public FilterModel
  {
public:
  Simple3T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(15, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(10, 10) = _Q(11, 11) = _Q(12, 12) = qs;

    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = _Q(13, 13) = _Q(14, 14) = ql;
  }

  virtual ~Simple3T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor3T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2, vec3_t& m3,
                            vec3_t& l3);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  };

#endif  // FILTER_MODEL_H_
