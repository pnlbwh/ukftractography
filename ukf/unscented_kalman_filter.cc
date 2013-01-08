/**
 * \file unscented_kalman_filter.cc
 * \brief implementation of unscented_kalman_filter.h
*/

#include "unscented_kalman_filter.h"
#include "filter_model.h"
#include "LU_solver.h"

#include <iostream>
#include <cassert>

#include "QuadProg++_vnl.h"

#include <limits>

#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_cholesky.h>

using namespace QuadProgPP;
using namespace LU_Solver;

UnscentedKalmanFilter::UnscentedKalmanFilter(FilterModel *filter_model)
  : _filter_model(filter_model), _k(0.01)
{
  int dim = _filter_model->state_dim();
  _scale = sqrt(dim + _k);

  _w.push_back(_k / (dim + _k));
  _w.insert(_w.end(), dim * 2, 0.5 / (dim + _k));

  // Create diagonal matrix.
  _w_.set_size(2 * dim + 1, 2 * dim + 1);
  _w_.fill(0); // clean memory left overs

  for (int i = 0; i < 2 * dim + 1; ++i) {
    _w_(i, i) = _w[i];
  }

  assert(static_cast<int>((_filter_model->Q()).rows()) == dim &&
         static_cast<int>((_filter_model->Q()).cols()) == dim);

  int signal_dim = _filter_model->signal_dim();
  assert(static_cast<int>((_filter_model->R()).rows()) == signal_dim &&
         static_cast<int>((_filter_model->R()).cols()) == signal_dim);

  // Resizing of temporary storage.
  // NOTE: filling all temp matrices with zeroes, to make sure they are really empty.
  dim_dimext.set_size(dim, 2 * dim + 1);
  dim_dimext.fill(0);
  signaldim_dimext.set_size(signal_dim, 2 * dim + 1);
  signaldim_dimext.fill(0);
  dim_dim.set_size(dim, dim);
  dim_dim.fill(0);
  dim_signaldim.set_size(dim, signal_dim);
  dim_signaldim.fill(0);
  signaldim_dim.set_size(signal_dim, dim);
  signaldim_dim.fill(0);
  X.set_size(dim, 2 * dim + 1);
  X.fill(0);
  X_.set_size(dim, 2 * dim + 1);
  X_.fill(0);
  Y.set_size(signal_dim, 2 * dim + 1);
  Y.fill(0);
  Pxy.set_size(dim, signal_dim);
  Pxy.fill(0);
  Pyy.set_size(signal_dim, signal_dim);
  Pyy.fill(0);
  K.set_size(signal_dim, dim);
  K.fill(0);

  X_hat.set_size(dim);
  Y_hat.set_size(signal_dim);
}

void UnscentedKalmanFilter::SigmaPoints(const State& x,
                                        const vnl_matrix<double>& p,
                                        vnl_matrix<double>& x_spread)
{

  int dim = _filter_model->state_dim();

  assert(x_spread.rows() == static_cast<unsigned int>(dim) &&
         x_spread.cols() == static_cast<unsigned int>(2 * dim + 1));
  assert(x.size() == static_cast<unsigned int>(dim));
  assert(p.rows() == static_cast<unsigned int>(dim) &&
         p.cols() == static_cast<unsigned int>(dim));

  State tmp_x = x;
  vnl_vector_ref<double> x_vnl(tmp_x.size(), &tmp_x.front());


  vnl_matrix<double> M(dim, dim);
  vnl_cholesky cholesky_decomp(p, vnl_cholesky::quiet);
  M = cholesky_decomp.lower_triangle();

  M = M * _scale;

  // Horizontally stack x to the X_tmp matrix.
  vnl_matrix<double> X_tmp(dim, dim); // CB: X changed to X_tmp to avoid notation confusion with member var X

  X.fill(0);
  for (size_t i = 0; i < X_tmp.cols(); ++i) {
    X_tmp.set_column(i, x_vnl);
  }

  // Create dim x (2 * dim + 1) matrix (x, x + m, x - m).
  x_spread.set_column(0, x_vnl);
  x_spread.update(X_tmp + M, 0, 1);
  x_spread.update(X_tmp - M, 0, dim + 1);

}

// vector version
void UnscentedKalmanFilter::Constrain(vnl_vector<double>& x, const vnl_matrix<double>& W)
{

  int dim_state = x.size();
  double error;

  // The equality constraints are just dummy variables. The solve_quadprog function has been changed
  // to ignore equality constraints.
  vnl_matrix<double> CE(dim_state, 1);
  CE.fill(0);
  vnl_vector<double> ce0(1);
  ce0.fill(0);

  vnl_matrix<double> W_tmp = W; 	 // W will be changed in solve_quadprog
  W_tmp = (W_tmp + W_tmp.transpose()) / 2; // ensure symmetry

  if (violatesContraints(x)) {
    vnl_vector<double> g0 = -1.0 * (W_tmp.transpose()) * x;
    vnl_vector<double> d  = _filter_model->d(); // the inequality constraints
    vnl_matrix<double> D  = _filter_model->D(); // -- " --
    error = solve_quadprog(W_tmp, g0, CE, ce0, D, d, x);
    if (error > 0.01) { // error usually much smaller than that, if solve_quadprog fails it returns inf
      std::cout << "Error too big while constraining the state. It was: " << error << std::endl;
      exit(1);
    }
  }
}

// matrix version
void UnscentedKalmanFilter::Constrain(vnl_matrix<double>& X, const vnl_matrix<double>& W)
{

  for (unsigned int i = 0; i < X.cols(); ++i) {
    vnl_vector<double> x = X.get_column(i);
    Constrain(x, W);
    X.set_column(i, &x[0]);
  }
}

bool UnscentedKalmanFilter::violatesContraints(vnl_vector<double> & x)
{
  vnl_vector<double> d_test = (-1.0) * ((_filter_model->D().transpose()) * x); 	// -D'*x
  for (unsigned int i = 0; i < d_test.size(); ++i) {						// if any(-D'*x > d) constraint is broken
    if (d_test[i]  > (_filter_model->d())[i]) {
      return true;
    }
  }
  return false;
}

void UnscentedKalmanFilter::Filter(const State& x,
                                   const vnl_matrix<double>& p,
                                   const std::vector<double>& z,
                                   State& x_new,
                                   vnl_matrix<double>& p_new,
                                   double& dNormMSE )
{

  assert(static_cast<int>(x.size()) == _filter_model->state_dim());
  assert(static_cast<int>(z.size()) == _filter_model->signal_dim());
  assert(static_cast<int>(x_new.size()) == _filter_model->state_dim());
  assert(static_cast<int>(p.rows()) == _filter_model->state_dim() &&
         static_cast<int>(p.cols()) == _filter_model->state_dim());
  assert(static_cast<int>(p_new.rows()) == _filter_model->state_dim() &&
         static_cast<int>(p_new.cols()) == _filter_model->state_dim());
  assert(_filter_model);
  assert(static_cast<int>(_w.size()) == 2 * _filter_model->state_dim() + 1);

  const vnl_matrix<double> Q = _filter_model->Q();
  const vnl_matrix<double> R = _filter_model->R();

  // Create sigma points.
  SigmaPoints(x, p, X); // doesnt change p, its const

  if (_filter_model->isConstrained()) {
    //vnl_matrix<double> p_tmp = p; // will be changed in QuadProg
    Constrain(X, p);
  }

  _filter_model->F(X); // slightly negative fw is fixed here

  vnl_vector_ref<double> _w_vnl(_w.size(), &_w.front());
  vnl_vector_ref<double> x_new_vnl(x_new.size(), &x_new.front());

  // copy step needed because z is a const variable and can't be referenced
  std::vector<double> z_tmp = z;
  vnl_vector_ref<double> z_vnl(z_tmp.size(), &z_tmp.front());

  X_hat = X * _w_vnl;
  for (size_t i = 0; i < dim_dimext.cols(); ++i) {
    dim_dimext.set_column(i, X_hat);
  }
  X_ = X - dim_dimext;

  p_new = X_ * _w_ * X_.transpose() + Q;

  _filter_model->H(X, Y);

  Y_hat = Y * _w_vnl;
  for (size_t i = 0; i < signaldim_dimext.cols(); ++i) {
    signaldim_dimext.set_column(i, Y_hat);
  }

  Y -= signaldim_dimext;
  vnl_matrix<double>& Y_ = Y;

  Pyy = Y_ * _w_ * Y_.transpose() + R;

  // Predict cross-correlation between state and observation.
  Pxy = X_ * _w_ * Y_.transpose();

  // Kalman gain K, estimate state/observation, compute covariance.
  // Solve K = Pyy \ Pxy'
  signaldim_dim = Pxy.transpose();

  vnl_matrix<double> LU = Pyy;
  LUdecmpDoolittle(&LU(0, 0), LU.rows());
  LUsolveDoolittle(&LU(0, 0), &signaldim_dim(0, 0), signaldim_dim.rows(), signaldim_dim.cols());

  K = signaldim_dim;

  dNormMSE = ((z_vnl - Y_hat).squared_magnitude()) / (z_vnl.squared_magnitude());

  p_new = p_new - K.transpose() * Pyy * K;
  Y_hat = z_vnl - Y_hat; // z is the real signal

  // NOTE: 	x_new_vnl = ... wont work for the reference type because the operator= is private.
  // 		Instead the copy in function of the vnl_vector base class is used.
  x_new_vnl.copy_in( &((K.transpose() * Y_hat + X_hat)[0]) );

  if (_filter_model->isConstrained()) {
    Constrain(x_new_vnl, p_new);
  }

}


