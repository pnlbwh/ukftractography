/**
 * \file unscented_kalman_filter.cc
 * \brief implementation of unscented_kalman_filter.h
*/

#include "unscented_kalman_filter.h"
#include "filter_model.h"

#include <iostream>
#include <cassert>


#include <limits>
#include <algorithm>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_cholesky.h>

#include "QuadProg++_Eigen.h"
using namespace QuadProgPP;
//#include "LU_solver.h"
//using namespace LU_Solver;

// HACK:  A quick conversion code for testing conversions between libraries.
namespace
{
inline Eigen::MatrixXd ToMatrixXd( const ukfMatrixType & in)
{
  Eigen::MatrixXd out;
  out.resize(in.rows(),in.cols());
  for(size_t r=0; r< in.rows(); ++r)
    {
    for(size_t c=0; c< in.cols(); ++c)
      {
      out(r,c)=in[r][c];
      }
    }
  return out;
}
inline ukfMatrixType ToVNL( const Eigen::MatrixXd & in)
{
  ukfMatrixType out;
  out.set_size(in.rows(),in.cols());
  for(int r=0; r< in.rows(); ++r)
    {
    for(int c=0; c< in.cols(); ++c)
      {
      out[r][c]=in(r,c);
      }
    }
  return out;
}

template <typename TIn, typename TOut>
TOut ConvertVector(const TIn &in)
{
  TOut out(in.size());
  //std::copy(in.begin(),in.end(),out.begin());
  for (unsigned int i =0 ; i < in.size(); ++i )
  {
    out[i] = in[i];
  }
  return out;
}

}

UnscentedKalmanFilter::UnscentedKalmanFilter(FilterModel *filter_model)
  : m_FilterModel(filter_model), m_SigmaPointSpread(0.01)
{
  const int dim = m_FilterModel->state_dim();

  m_Scale = sqrt(dim + m_SigmaPointSpread);

  m_Weights.push_back(m_SigmaPointSpread / (dim + m_SigmaPointSpread) );
  m_Weights.insert(m_Weights.end(), dim * 2, 0.5 / (dim + m_SigmaPointSpread) );

  // HANS NOTE: I think this would be the equivalent for Eigen
  // m_Weights.resize((2 * dim) + 1);
  // m_Weights(0) = m_SigmaPointSpread / (dim + m_SigmaPointSpread);
  // for(unsigned i = 1; i < (2 * dim) + 1; ++i)
  //   {
  //   m_Weights(i) = 0.5 / (dim + m_SigmaPointSpread);
  //   }

  // Create diagonal matrix.
  m_WeightsRepeated.resize(2 * dim + 1,2 * dim + 1);
  m_WeightsRepeated.setConstant(0.0);
  for( int i = 0; i < 2 * dim + 1; ++i )
    {
    m_WeightsRepeated(i, i) = m_Weights[i];
    }

  assert(static_cast<int>( (m_FilterModel->Q() ).rows() ) == dim &&
         static_cast<int>( (m_FilterModel->Q() ).cols() ) == dim);
}

void UnscentedKalmanFilter::SigmaPoints(const State& x,
                                        const Eigen::MatrixXd& p,
                                        Eigen::MatrixXd& x_spread)
{
  const unsigned int dim = m_FilterModel->state_dim();

  assert(x_spread.rows() == dim &&
         x_spread.cols() == static_cast<unsigned int>(2 * dim + 1) );
  assert(x.size() == dim );
  assert(p.rows() == dim &&
         p.cols() == dim );

  const Eigen::VectorXd x_vnl = ConvertVector<State,Eigen::VectorXd >(x);

  // Horizontally stack x to the X_tmp matrix.
  Eigen::MatrixXd X_tmp(dim, dim); // CB: X changed to X_tmp to avoid notation confusion with member var X
  for( unsigned int c = 0; c < dim; ++c )
    {
    X_tmp.col(c) = x_vnl;
    }
  Eigen::LLT<Eigen::MatrixXd> lltOfA(p); // compute the Cholesky decomposition of A
  Eigen::MatrixXd NewM = ( lltOfA.matrixL() ); // retrieve factor L  in the decomposition
  NewM *= m_Scale;

  // Create dim x (2 * dim + 1) matrix (x, x + m, x - m).
  x_spread.col(0) = x_vnl;
  x_spread.block(0,    1,dim,dim) = X_tmp + NewM;
  x_spread.block(0,dim+1,dim,dim) = X_tmp - NewM;
}

// vector version
void UnscentedKalmanFilter::Constrain(Eigen::VectorXd& x, const Eigen::MatrixXd& W)
{
  if( violatesContraints(x) )
    {
    const int    dim_state = x.size();
    // The equality constraints are just dummy variables. The solve_quadprog function has been changed
    // to ignore equality constraints.
    Eigen::MatrixXd CE(dim_state, 1);
    CE.setConstant(0.0);
    Eigen::VectorXd ce0(1);
    ce0.setConstant(0.0);

    const Eigen::MatrixXd WTranspose = W.transpose();
    //TODO:  review solve_quadprog to determine if any of these are const parameters.
    Eigen::MatrixXd W_tmp = (W + WTranspose) * 0.5;

    Eigen::VectorXd g0 = -1.0 * (W_tmp.transpose() ) * x;
    const Eigen::VectorXd d  = ConvertVector< ukfVectorType, Eigen::VectorXd>( m_FilterModel->d() ); // the inequality constraints
    const Eigen::MatrixXd D  = ToMatrixXd( m_FilterModel->D() ); // -- " --
    const double error = solve_quadprog(W_tmp, g0, CE, ce0, D, d, x);
    if( error > 0.01 )   // error usually much smaller than that, if solve_quadprog fails it returns inf
      {
      exit(1);
      }
    }
}

// matrix version
void UnscentedKalmanFilter::Constrain(Eigen::MatrixXd& localX, const Eigen::MatrixXd& localW)
{
  const unsigned int numCols = localX.cols();

  for( unsigned int i = 0; i < numCols; ++i )
    {
    Eigen::VectorXd x = localX.col(i);
    Constrain(x, localW);
    localX.col(i) = x;
    }
}

bool UnscentedKalmanFilter::violatesContraints(Eigen::VectorXd & x)
{
  const Eigen::VectorXd d_test = (-1.0) * ( ToMatrixXd( m_FilterModel->D().transpose() ) * x); // -D'*x
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
                                   const ukfMatrixType& p_VNL,
                                   const std::vector<double>& z,
                                   State& x_new,
                                   ukfMatrixType& p_new,
                                   double& dNormMSE )
{
  const Eigen::MatrixXd p = ToMatrixXd( p_VNL );
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
  Eigen::MatrixXd X(dim, 2 * dim + 1);
  X.setConstant(0.0);
  // Create sigma points.
  SigmaPoints(x, p, X); // doesnt change p, its const

  if( localConstFilterModel->isConstrained() )
    {
    // Eigen::MatrixXd p_tmp = p; // will be changed in QuadProg
    Constrain(X, p);
    }

  {
  ukfMatrixType X_toFix = ToVNL(X); //HACK
  localConstFilterModel->F(X_toFix); // slightly negative fw is fixed here
  X = ToMatrixXd( X_toFix );
  }
  // copy step needed because z is a const variable and can't be referenced
  Eigen::VectorXd z_vnl(z.size());
  for( size_t i = 0; i < z.size(); ++i)
  {
    z_vnl[i]=z[i];
  }

  /** Used for the estimation of the new state */
  Eigen::MatrixXd dim_dimext(dim, 2 * dim + 1);
  dim_dimext.setConstant(0.0);
  Eigen::VectorXd _w_vnl = ConvertVector<std::vector<double>, Eigen::VectorXd >(this->m_Weights);
  const Eigen::VectorXd X_hat = X * _w_vnl;
  for( unsigned int i = 0; i < dim_dimext.cols(); ++i )
    {
    dim_dimext.col(i) = X_hat;
    }

  const Eigen::MatrixXd & Q = ToMatrixXd( localConstFilterModel->Q() );

  const Eigen::MatrixXd &X_ = X - dim_dimext;
  // Use const reference to avoid copying
  Eigen::MatrixXd p_new_Eigen = X_ * m_WeightsRepeated * X_.transpose() + Q;

  /** The signal */
  Eigen::MatrixXd Y(signal_dim, 2 * dim + 1);
  Y.setConstant(0.0);
    { // HACK
    ukfMatrixType X_VNL = ToVNL( X );
    ukfMatrixType Y_VNL = ToVNL( Y );
    localConstFilterModel->H(X_VNL, Y_VNL);
    X = ToMatrixXd( X_VNL );
    Y = ToMatrixXd( Y_VNL );
    }

  /** Used for the estimation of the signal */

  Eigen::MatrixXd     signaldim_dimext(signal_dim, 2 * dim + 1);
  signaldim_dimext.setConstant(0.0);
  const Eigen::VectorXd Y_hat = Y * _w_vnl;
  for( unsigned int i = 0; i < signaldim_dimext.cols(); ++i )
    {
    signaldim_dimext.col(i) = Y_hat;
    }

  Y -= signaldim_dimext;
  const Eigen::MatrixXd Y_ = Y;
  const Eigen::MatrixXd WeightsRepeated_Y_Transpose =m_WeightsRepeated*Y_.transpose();

  const Eigen::MatrixXd R = ToMatrixXd( localConstFilterModel->R() );

  /** Covariance of the signal */
  const Eigen::MatrixXd Pyy = Y_ * WeightsRepeated_Y_Transpose + R;

  // Predict cross-correlation between state and observation.
  /** Covariance matrix state/signal */
  const Eigen::MatrixXd Pxy = X_ * WeightsRepeated_Y_Transpose;

  // Kalman gain KalmanGainMatrix, estimate state/observation, compute covariance.
  // Solve KalmanGainMatrix = Pyy \ Pxy'
  // Solve Ax = b. Result stored in x. Matlab: x = A \ b.
  //x = A.ldlt().solve(b));
  const Eigen::MatrixXd KalmanGainMatrix = Pyy.ldlt().solve(Pxy.transpose());

  dNormMSE = ( (z_vnl - Y_hat).squaredNorm() ) / (z_vnl.squaredNorm() );

  p_new_Eigen = p_new_Eigen - KalmanGainMatrix.transpose() * Pyy * KalmanGainMatrix;
  const Eigen::VectorXd Y_hat2 = z_vnl - Y_hat; // z is the real signal

  // HANS NOTE -- the original code (through the vector_ref)
  // over_wrote x_new with the following expression's result, then
  // optionally called Constrain on it. So this code doesn't bother
  // copying the old value of x_new into x_new_vnl before doing this
  // computation, it just copies the resut in afterwards.
  Eigen::VectorXd x_new_Eigen = KalmanGainMatrix.transpose() * Y_hat2 + X_hat;
  if( localConstFilterModel->isConstrained() )
    {
    Constrain(x_new_Eigen, p_new_Eigen);
    }
  x_new = ConvertVector<Eigen::VectorXd,State>(x_new_Eigen);
  p_new = ToVNL( p_new_Eigen );
}
