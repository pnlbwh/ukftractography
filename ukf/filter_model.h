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
struct FilterModel
  {

  /** Constructor */
  FilterModel(int state_dim, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
    : _state_dim(state_dim), _rs(rs), _signal_dim(0), _signal_data(NULL), weights_on_tensors_(weights_on_tensors),
    _constrained(constrained)
  {

    _Q.set_size(_state_dim, _state_dim);
    _Q.fill(0); // necessary because otherwise there is memory left overs in the matrix

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
  virtual void F(vnl_matrix<double>& X) const = 0;

  /** observation, i.e. signal reconstruction */
  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const = 0;

  // Functions that convert the state into a tensor representation that can be
  // used for tractography (meaning the main direction and all the eigenvalues/
  // lambdas).

  /** Extracts principal diffusion direction and eigen values from the state for the 1T cases */
  virtual void State2Tensor(const State &, vec_t &, vec_t & )
  {
    assert(!"Not implemented");
  }

  /** Extracts principal diffusion direction and eigen values from the state for the 2T cases */
  virtual void State2Tensor(const State &, const vec_t &, vec_t &,
                            vec_t &, vec_t &, vec_t & )
  {
    assert(!"Not implemented");
  }

  /** Extracts principal diffusion direction and eigen values from the state for the 3T cases */
  virtual void State2Tensor(const State &, const vec_t &, vec_t &,
                            vec_t &, vec_t &, vec_t &, vec_t &,
                            vec_t & )
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

    _R.set_size(_signal_dim, _signal_dim);
    _R.fill(0); // necessary because otherwise there is memory leftovers in the matrix
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
  const vnl_matrix<double> & Q() const
  {
    return _Q;
  }

  /** The noise in the signal reconstruction used by the UKF. */
  const vnl_matrix<double> & R() const
  {
    return _R;
  }

  // The next two functions are only used for the constrained case

  /** The inequality constraint matrix for the constrained UKF */
  const vnl_matrix<double> & D() const
  {
    return _D;
  }

  /** The inequality constraint right hand side for the constrained UKF */
  const vnl_vector<double> & d() const
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

  /** Checks if d is smaller than a small negative threshold. If yes an error is returned. Otherwise d is rounded to 0.0
    */
  double CheckZero(const double & d) const;

  /** The dimension of the state */
  const int _state_dim;

  /** The constant signal noise for each signal component */
  double _rs;

  /** The dimension of the signal */
  int _signal_dim;

  /** Pointer to the diffusion signal data */
  ISignalData *_signal_data;

  /** Process noise */
  vnl_matrix<double> _Q;

  /** Signal reconstruction noise */
  vnl_matrix<double> _R;

  /** Inequality constraint matrix, only used for constrained UKF */
  vnl_matrix<double> _D;

  /** Inequality right hand side, only used for constrained UKF */
  vnl_vector<double> _d;

  /** The weights of each tensor. (Most commonly equal all are equal) */
  const std::vector<double> weights_on_tensors_;

  /** Are we using the constrained filter */
  bool _constrained;
  };

/**
 * \struct Full1T_FW
 * \brief Full 1-Tensor model with free water
 *
 * Model describing 1-tensor tractography with the full tensor representation (3 angles, 3 eigenvalues)
 * and free water estimation.
*/
struct Full1T_FW : public FilterModel
  {
  Full1T_FW(double qs, double ql, double qw, double rs, const std::vector<double>& weights_on_tensors, bool constrained,
            const double diff_fw)
    : FilterModel(7, rs, weights_on_tensors, constrained), _lambda_min(100.0), _d_iso(diff_fw)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
    _Q(6, 6) = qw;

    _D.set_size(7, 5);
    _D.fill(0);

    _d.set_size(5);

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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, vec_t& m, vec_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  /** apparent diffusion coefficient of free water */
  const double _d_iso;

  };

/**
 * \struct Full1T
 * \brief Full 1-Tensor model without free water
 *
 * Model describing 1-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
struct Full1T : public FilterModel
  {
  Full1T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
  }

  virtual ~Full1T()
  {
  }

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, vec_t& m, vec_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;
  };

/**
 * \struct Full2T
 * \brief Full 2-Tensor model without free water
 *
 * Model describing 2-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
struct Full2T : public FilterModel
  {
  Full2T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
    : FilterModel(12, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
  }

  virtual ~Full2T()
  {
  }

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;
  };

/**
 * \struct Full2T_FW
 * \brief Full 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the full tensor representation (3 angles, 3 eigenvalues)
 * and free water estimation.
*/
struct Full2T_FW : public FilterModel
  {
  Full2T_FW(double qs, double ql, double qw, double rs, const std::vector<double>& weights_on_tensors, bool constrained,
            const double diff_fw)
    : FilterModel(13, rs, weights_on_tensors, constrained), _lambda_min(100.0), _d_iso(diff_fw)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
    _Q(12, 12) = qw; // noise for weights

    _D.set_size(13, 8);
    _D.fill(0);

    _d.set_size(8);

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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  /** apparent diffusion coefficient of free water */
  const double _d_iso;

  };

/**
 * \struct Full3T
 * \brief Full 3-Tensor model with free water
 *
 * Model describing 3-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
struct Full3T : public FilterModel
  {
  Full3T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2, vec_t& m3,
                            vec_t& l3);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  };

/**
 * \struct Simple1T_FW
 * \brief Simple 1-Tensor model with free water
 *
 * Model describing 1-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
struct Simple1T_FW : public FilterModel
  {
  Simple1T_FW(double qs, double ql, double qw, double rs, const std::vector<double>& weights_on_tensors,
              bool constrained, const double diff_fw)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0), _d_iso(diff_fw)
  {

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = ql;
    _Q(5, 5) = qw; // noise for weights

    _D.set_size(6, 4);
    _D.fill(0);

    _d.set_size(4);

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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, vec_t& m, vec_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  /** apparent diffusion coefficient of free water */
  const double _d_iso;

  };

/**
 * \struct Simple1T
 * \brief Simple 1-Tensor model
 *
 * Model describing 1-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
*/
struct Simple1T : public FilterModel
  {
  Simple1T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
    : FilterModel(5, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = ql;
  }

  virtual ~Simple1T()
  {
  }

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, vec_t& m, vec_t& l);

  const double _lambda_min;
  };

/**
 * \struct Simple2T
 * \brief Simple 2-Tensor model
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal).
*/
struct Simple2T : public FilterModel
  {
  Simple2T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
    : FilterModel(10, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
  }

  virtual ~Simple2T()
  {
  }

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  };

/**
 * \struct Simple2T_FW
 * \brief Simple 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
struct Simple2T_FW : public FilterModel
  {
  Simple2T_FW(double qs, double ql, double qw, double rs, const std::vector<double>& weights_on_tensors,
              bool constrained, const double diff_fw)
    : FilterModel(11, rs, weights_on_tensors, constrained), _lambda_min(100.0), _d_iso(diff_fw)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
    _Q(10, 10) = qw; // noise for weights

    _D.set_size(11, 6);
    _D.fill(0);

    _d.set_size(6);

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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  /** apparent diffusion coefficient of free water */
  const double _d_iso;
  };

/**
 * \struct Simple3T
 * \brief Simple 3-Tensor model
 *
 * Model describing 3-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
*/
struct Simple3T : public FilterModel
  {
  Simple3T(double qs, double ql, double rs, const std::vector<double>& weights_on_tensors, bool constrained)
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

  virtual void F(vnl_matrix<double>& X) const;

  virtual void H(const  vnl_matrix<double>& X, vnl_matrix<double>& Y) const;

  virtual void State2Tensor(const State& x, const vec_t& old_m, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2, vec_t& m3,
                            vec_t& l3);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const double _lambda_min;

  };

#endif  // FILTER_MODEL_H_
