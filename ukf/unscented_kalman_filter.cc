/**
 * \file unscented_kalman_filter.cc
 * \brief implementation of unscented_kalman_filter.h
*/

#include "unscented_kalman_filter.h"
#include "filter_model.h"

#include <iostream>
#include <cassert>
#include <stdexcept>

#include <limits>
#include <algorithm>

#include "QuadProg++_Eigen.h"
using namespace QuadProgPP;
//#include "LU_solver.h"
//using namespace LU_Solver;

UnscentedKalmanFilter::UnscentedKalmanFilter(FilterModel *filter_model)
  : m_FilterModel(filter_model), m_SigmaPointSpread(0.01)
{
  const unsigned int dim = m_FilterModel->state_dim();

  m_Scale = sqrt(dim + m_SigmaPointSpread);
  m_Weights.resize((2 * dim) + 1);
  m_Weights(0) = m_SigmaPointSpread / (dim + m_SigmaPointSpread);
  // Create diagonal matrix.
  m_WeightsRepeated.resize(2 * dim + 1,2 * dim + 1);
  m_WeightsRepeated.setConstant(ukfZero);
  m_WeightsRepeated(0, 0) = m_Weights[0];
  for(unsigned int i = 1; i < (2 * dim) + 1; ++i)
     {
     m_Weights(i) = ukfHalf / (dim + m_SigmaPointSpread);
     m_WeightsRepeated(i, i) = m_Weights[i];
     }

  assert(static_cast<unsigned int>( (m_FilterModel->Q() ).rows() ) == dim &&
         static_cast<unsigned int>( (m_FilterModel->Q() ).cols() ) == dim);

  //DUMMY VARIABLES TAHT ARE ALWAYS ZERO and not used.
  m_DummyZeroCE.resize(dim, 1);
  m_DummyZeroCE.setConstant(ukfZero);
  m_DummyZeroce0.resize(1);
  m_DummyZeroce0.setConstant(ukfZero);
}

void UnscentedKalmanFilter::SigmaPoints(const State& x,
                                        const ukfMatrixType& p,
                                        ukfMatrixType& x_spread)
{
  const unsigned int dim = m_FilterModel->state_dim();

  assert(x_spread.rows() == dim &&
         x_spread.cols() == static_cast<unsigned int>(2 * dim + 1) );
  assert(x.size() == dim );
  assert(p.rows() == dim &&
         p.cols() == dim );

  // Horizontally stack x to the X_tmp matrix.
  ukfMatrixType X_tmp(dim, dim); // CB: X changed to X_tmp to avoid notation confusion with member var X
  for( unsigned int c = 0; c < dim; ++c )
    {
    X_tmp.col(c) = x;
    }
  Eigen::LLT<ukfMatrixType> lltOfA(p); // compute the Cholesky decomposition of A
  ukfMatrixType NewM = ( lltOfA.matrixL() ); // retrieve factor L  in the decomposition
  NewM *= m_Scale;
  // Create dim x (2 * dim + 1) matrix (x, x + m, x - m).
  x_spread.col(0) = x;
  x_spread.block(0,    1,dim,dim) = X_tmp + NewM;
  x_spread.block(0,dim+1,dim,dim) = X_tmp - NewM;
}

// vector version
void UnscentedKalmanFilter::Constrain(ukfVectorType& x, const ukfMatrixType& W)
{
  if( violatesContraints(x) )
    {
    const ukfMatrixType WTranspose = W.transpose();
    ukfMatrixType W_tmp = (W + WTranspose) * ukfHalf;
    ukfVectorType g0 = -ukfOne * (W_tmp.transpose() ) * x;
    const ukfVectorType d  = m_FilterModel->d(); // the inequality constraints
    const ukfMatrixType D  = m_FilterModel->D(); // -- " --
    // The equality constraints are just dummy variables. The solve_quadprog function has been changed
    // to ignore equality constraints.
    const ukfPrecisionType error = solve_quadprog(W_tmp, g0, m_DummyZeroCE, m_DummyZeroce0, D, d, x);
    if( error > 0.01 )   // error usually much smaller than that, if solve_quadprog fails it returns inf
      {
      throw std::logic_error("solve_quadprog error exceeds threshold 0.01!");
      }
    }
}

// matrix version
void UnscentedKalmanFilter::Constrain(ukfMatrixType& localX, const ukfMatrixType& localW)
{
  const unsigned int numCols = static_cast<unsigned int>(localX.cols());

  for( unsigned int i = 0; i < numCols; ++i )
    {
    ukfVectorType x = localX.col(i);
    Constrain(x, localW);
    localX.col(i) = x;
    }
}

bool UnscentedKalmanFilter::violatesContraints(ukfVectorType & x)
{
  const ukfVectorType d_test = (-ukfOne) * ( m_FilterModel->D().transpose() ) * x; // -D'*x
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
                                   const ukfVectorType& z,
                                   State& x_new,
                                   ukfMatrixType& p_new,
                                   ukfPrecisionType& dNormMSE )
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

  const int dim = localConstFilterModel->state_dim();
  const int signal_dim = localConstFilterModel->signal_dim();

  /** The state spread out according to the unscented transform */
  ukfMatrixType X(dim, 2 * dim + 1);
  X.setConstant(ukfZero);
  // Create sigma points.
  SigmaPoints(x, p, X); // doesnt change p, its const
  // std::cout <<"\n state:"<<x;
  if( localConstFilterModel->isConstrained() )
    {
    // ukfMatrixType p_tmp = p; // will be changed in QuadProg
    Constrain(X, p);
    }

  {
  localConstFilterModel->F(X); // slightly negative fw is fixed here
  }
  // copy step needed because z is a const variable and can't be referenced
  ukfVectorType z_Eigen(z.size());
  for( unsigned int i = 0; i < z.size(); ++i)
  {
    z_Eigen[i]=z[i];
  }
  // std::cout<<"recorded signal: "<<z_Eigen;

  /** Used for the estimation of the new state */
  ukfMatrixType dim_dimext(dim, 2 * dim + 1);
  dim_dimext.setConstant(ukfZero);
  // std::cout << "\n X:"<<X <<"\n Weigths" << this->m_Weights;
  const ukfVectorType X_hat = X * this->m_Weights;
  for( unsigned int i = 0; i < dim_dimext.cols(); ++i )
    {
    dim_dimext.col(i) = X_hat;
    }

  const ukfMatrixType & Q = localConstFilterModel->Q();

  const ukfMatrixType &X_ = X - dim_dimext;
  // Use const reference to avoid copying
  p_new = X_ * m_WeightsRepeated * X_.transpose() + Q;
  // std::cout<<"\n intitial P:"<<p_new;
  /** The signal */
  const ukfMatrixType Yk = p_new.inverse();
  const ukfMatrixType yhat = Yk * X_hat;
  // std::cout << "\n P: " << p_new << "\n yhat: " << yhat <<"\n Q: "<<Q;
  // exit(1);
  ukfMatrixType Z(signal_dim, 2 * dim + 1);
  Z.setConstant(ukfZero);
  localConstFilterModel->H(X, Z);
  /** Used for the estimation of the signal */

  ukfMatrixType     signaldim_dimext(signal_dim, 2 * dim + 1);
  signaldim_dimext.setConstant(ukfZero);
  const ukfVectorType Z_hat = Z * this->m_Weights;
  for( unsigned int i = 0; i < signaldim_dimext.cols(); ++i )
    {
    signaldim_dimext.col(i) = Z_hat;
    }

  Z -= signaldim_dimext;
  const ukfMatrixType Z_ = Z;
  //const ukfMatrixType WeightsRepeated_Y_Transpose =m_WeightsRepeated*Y_.transpose();

  const ukfMatrixType temp = localConstFilterModel->R();
  ukfPrecisionType R = 1/temp(0,0);

  /** Covariance of the signal */
  //const ukfMatrixType Pyy = Y_ * WeightsRepeated_Y_Transpose + R;

  // Predict cross-correlation between state and observation.
  /** Covariance matrix state/signal */
  const ukfMatrixType Pxz = X_ * m_WeightsRepeated * Z_.transpose();

  // Kalman gain KalmanGainMatrix, estimate state/observation, compute covariance.
  // Solve KalmanGainMatrix = Pyy \ Pxy'
  // Solve Ax = b. Result stored in x. Matlab: x = A \ b.
  //x = A.ldlt().solve(b));
  const ukfMatrixType Ht = Yk * Pxz;

  dNormMSE = ( (z_Eigen - Z_hat).squaredNorm() ) / (z_Eigen.squaredNorm() );

  const ukfMatrixType I = Ht * R * Ht.transpose();
  const ukfMatrixType i = Ht * R * ((z_Eigen - Z_hat) + (Pxz.transpose() * yhat));
  // std::cout <<"\n I:"<<I<<"\n i :"<<i;
  p_new = (Yk+I).inverse();
  // std::cout <<"\n p_new:"<<p_new;
  ukfVectorType x_new_Eigen = p_new * (i + yhat);
  if( localConstFilterModel->isConstrained() )
    {
    Constrain(x_new_Eigen, Yk+I);
    }
  x_new = x_new_Eigen;
}
