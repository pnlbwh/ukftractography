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
  : m_FilterModel(filter_model), m_SigmaPointSpread(0.01)
{
  const int dim = m_FilterModel->state_dim();

  m_Scale = sqrt(dim + m_SigmaPointSpread);

  m_Weights.push_back(m_SigmaPointSpread / (dim + m_SigmaPointSpread) );
  m_Weights.insert(m_Weights.end(), dim * 2, 0.5 / (dim + m_SigmaPointSpread) );

  // Create diagonal matrix.
  m_WeightsRepeated.set_size(2 * dim + 1, 2 * dim + 1);
  m_WeightsRepeated.fill(0); // clean memory left overs
  for( int i = 0; i < 2 * dim + 1; ++i )
    {
    m_WeightsRepeated(i, i) = m_Weights[i];
    }

  assert(static_cast<int>( (m_FilterModel->Q() ).rows() ) == dim &&
         static_cast<int>( (m_FilterModel->Q() ).cols() ) == dim);

  int signal_dim = m_FilterModel->signal_dim();
  assert(static_cast<int>( (m_FilterModel->R() ).rows() ) == signal_dim &&
         static_cast<int>( (m_FilterModel->R() ).cols() ) == signal_dim);

  // Resizing of temporary storage.
  // NOTE: filling all temp matrices with zeroes, to make sure they are really empty.
  m_KalmanGainMatrix.set_size(signal_dim, dim);
  m_KalmanGainMatrix.fill(0);

}

void UnscentedKalmanFilter::SigmaPoints(const State& x,
                                        const ukfMatrixType& p,
                                        ukfMatrixType& x_spread)
{
  const unsigned int dim = m_FilterModel->state_dim();

  assert(x_spread.rows() == static_cast<unsigned int>(dim) &&
         x_spread.cols() == static_cast<unsigned int>(2 * dim + 1) );
  assert(x.size() == static_cast<unsigned int>(dim) );
  assert(p.rows() == static_cast<unsigned int>(dim) &&
         p.cols() == static_cast<unsigned int>(dim) );

  State                  tmp_x = x;
  vnl_vector_ref<double> x_vnl(tmp_x.size(), &tmp_x.front() );

  vnl_cholesky       cholesky_decomp(p, vnl_cholesky::quiet);
  const ukfMatrixType &M = cholesky_decomp.lower_triangle()*m_Scale;

  // Horizontally stack x to the X_tmp matrix.
  ukfMatrixType X_tmp(dim, dim); // CB: X changed to X_tmp to avoid notation confusion with member var X
  for( unsigned int i = 0; i < dim; ++i )
    {
    X_tmp.set_column(i, x_vnl);
    }

  // Create dim x (2 * dim + 1) matrix (x, x + m, x - m).
  x_spread.set_column(0, x_vnl);
  x_spread.update(X_tmp + M, 0, 1);
  x_spread.update(X_tmp - M, 0, dim + 1);

}

// vector version
void UnscentedKalmanFilter::Constrain(ukfVectorType& x, const ukfMatrixType& W)
{
  if( violatesContraints(x) )
    {
    const int    dim_state = x.size();

    // The equality constraints are just dummy variables. The solve_quadprog function has been changed
    // to ignore equality constraints.
    ukfMatrixType CE(dim_state, 1);
    CE.fill(0);
    ukfVectorType ce0(1);
    ce0.fill(0);

#if 0
    ukfMatrixType W_tmp = W;             // W will be changed in solve_quadprog
    W_tmp = (W_tmp + W_tmp.transpose() ) / 2; // ensure symmetry
#else
    ukfMatrixType W_tmp = (W + W.transpose() ) *0.5;
#endif

    ukfVectorType g0 = -1.0 * (W_tmp.transpose() ) * x;
    ukfVectorType d  = m_FilterModel->d(); // the inequality constraints
    ukfMatrixType D  = m_FilterModel->D(); // -- " --
    const double error = solve_quadprog(W_tmp, g0, CE, ce0, D, d, x);
    if( error > 0.01 )   // error usually much smaller than that, if solve_quadprog fails it returns inf
      {
      std::cout << "Error too big while constraining the state. It was: " << error << std::endl;
      exit(1);
      }
    }
}

// matrix version
void UnscentedKalmanFilter::Constrain(ukfMatrixType& localX, const ukfMatrixType& localW)
{
  const unsigned int numCols = localX.cols();

  for( unsigned int i = 0; i < numCols; ++i )
    {
    ukfVectorType x = localX.get_column(i);
    Constrain(x, localW);
    localX.set_column(i, &x[0]);
    }
}

bool UnscentedKalmanFilter::violatesContraints(ukfVectorType & x)
{
  const ukfVectorType & d_test = (-1.0) * ( (m_FilterModel->D().transpose() ) * x); // -D'*x
  for( unsigned int i = 0; i < d_test.size(); ++i )                              // if any(-D'*x > d) constraint is
                                                                                 // broken
    {
    if( d_test[i]  > (m_FilterModel->d() )[i] )
      {
      return true;
      }
    }
  return false;
}

void UnscentedKalmanFilter::Filter(const State& x,
                                   const ukfMatrixType& p,
                                   const std::vector<double>& z,
                                   State& x_new,
                                   ukfMatrixType& p_new,
                                   double& dNormMSE )
{
  // Force a const version of the m_FilterModel to be used to ensure that it is not modified.
  FilterModel const * const localConstFilterModel = m_FilterModel;

  assert(static_cast<int>(x.size() ) == localConstFilterModel->state_dim() );
  assert(static_cast<int>(z.size() ) == localConstFilterModel->signal_dim() );
  assert(static_cast<int>(x_new.size() ) == localConstFilterModel->state_dim() );
  assert(static_cast<int>(p.rows() ) == localConstFilterModel->state_dim() &&
         static_cast<int>(p.cols() ) == localConstFilterModel->state_dim() );
  assert(static_cast<int>(p_new.rows() ) == localConstFilterModel->state_dim() &&
         static_cast<int>(p_new.cols() ) == localConstFilterModel->state_dim() );
  assert(localConstFilterModel);
  assert(static_cast<int>(m_Weights.size() ) == 2 * localConstFilterModel->state_dim() + 1);

  int dim = m_FilterModel->state_dim();
  int signal_dim = m_FilterModel->signal_dim();

  ukfMatrixType     signaldim_dimext;
  signaldim_dimext.set_size(signal_dim, 2 * dim + 1);
  signaldim_dimext.fill(0);

  ukfMatrixType     signaldim_dim;
  signaldim_dim.set_size(signal_dim, dim);
  signaldim_dim.fill(0);

  ukfMatrixType dim_dimext;
  dim_dimext.set_size(dim, 2 * dim + 1);
  dim_dimext.fill(0);

  ukfMatrixType dim_dim;
  dim_dim.set_size(dim, dim);
  dim_dim.fill(0);

  // Temporary storage.
  ukfMatrixType X_;
  X_.set_size(dim, 2 * dim + 1);
  X_.fill(0);

  /** The state spread out according to the unscented transform */
  ukfMatrixType X;
  X.set_size(dim, 2 * dim + 1);
  X.fill(0);

  /** The signal */
  ukfMatrixType Y;
  Y.set_size(signal_dim, 2 * dim + 1);
  Y.fill(0);

  /** Covariance matrix state/signal */
  ukfMatrixType Pxy;
  Pxy.set_size(dim, signal_dim);
  Pxy.fill(0);

  /** Covariance of the signal */
  ukfMatrixType Pyy;
  Pyy.set_size(signal_dim, signal_dim);
  Pyy.fill(0);

  /** Used for the estimation of the new state */
  ukfVectorType X_hat;
  X_hat.set_size(dim);

  /** Used for the estimation of the signal */
  ukfVectorType Y_hat;
  Y_hat.set_size(signal_dim);



  // Create sigma points.
  SigmaPoints(x, p, X); // doesnt change p, its const

  if( localConstFilterModel->isConstrained() )
    {
    // ukfMatrixType p_tmp = p; // will be changed in QuadProg
    Constrain(X, p);
    }

  localConstFilterModel->F(X); // slightly negative fw is fixed here

  // copy step needed because z is a const variable and can't be referenced
  std::vector<double>    z_tmp = z;
  vnl_vector_ref<double> z_vnl(z_tmp.size(), &z_tmp.front() );

  vnl_vector_ref<double> _w_vnl(m_Weights.size(), &m_Weights.front() );
  X_hat = X * _w_vnl;
  for( size_t i = 0; i < dim_dimext.cols(); ++i )
    {
    dim_dimext.set_column(i, X_hat);
    }
  X_ = X - dim_dimext;

  // Use const reference to avoid copying
  const ukfMatrixType & Q = localConstFilterModel->Q();
  p_new = X_ * m_WeightsRepeated * X_.transpose() + Q;

  localConstFilterModel->H(X, Y);

  Y_hat = Y * _w_vnl;
  for( size_t i = 0; i < signaldim_dimext.cols(); ++i )
    {
    signaldim_dimext.set_column(i, Y_hat);
    }

  Y -= signaldim_dimext;
  const ukfMatrixType& Y_ = Y;
  const ukfMatrixType& WeightsRepeated_Y_Transpose =m_WeightsRepeated*Y_.transpose();

  const ukfMatrixType & R = localConstFilterModel->R();
  Pyy = Y_ * WeightsRepeated_Y_Transpose + R;

  // Predict cross-correlation between state and observation.
  Pxy = X_ * WeightsRepeated_Y_Transpose;

  // Kalman gain m_KalmanGainMatrix, estimate state/observation, compute covariance.
  // Solve m_KalmanGainMatrix = Pyy \ Pxy'
  signaldim_dim = Pxy.transpose();

  ukfMatrixType LU = Pyy;
  LUdecmpDoolittle(&LU(0, 0), LU.rows() );
  LUsolveDoolittle(&LU(0, 0), &signaldim_dim(0, 0), signaldim_dim.rows(), signaldim_dim.cols() );

  m_KalmanGainMatrix = signaldim_dim;

  dNormMSE = ( (z_vnl - Y_hat).squared_magnitude() ) / (z_vnl.squared_magnitude() );

  p_new = p_new - m_KalmanGainMatrix.transpose() * Pyy * m_KalmanGainMatrix;
  Y_hat = z_vnl - Y_hat; // z is the real signal

  vnl_vector_ref<double> x_new_vnl(x_new.size(), &x_new.front() );
  // NOTE:  x_new_vnl = ... wont work for the reference type because the operator= is private.
  //    Instead the copy in function of the vnl_vector base class is used.
  x_new_vnl.copy_in( &( (m_KalmanGainMatrix.transpose() * Y_hat + X_hat)[0]) );

  if( localConstFilterModel->isConstrained() )
    {
    Constrain(x_new_vnl, p_new);
    }

}
