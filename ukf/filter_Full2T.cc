#include "filter_Full2T.h"

// Functions for 2-tensor full model.
void Full2T::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  // Clamp lambdas.
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    X(3, i) = std::max(X(3, i), _lambda_min);
    X(4, i) = std::max(X(4, i), _lambda_min);
    X(5, i) = std::max(X(5, i), _lambda_min);
    X(9, i) = std::max(X(9, i), _lambda_min);
    X(10, i) = std::max(X(10, i), _lambda_min);
    X(11, i) = std::max(X(11, i), _lambda_min);
    }
}

void Full2T::H(const  ukfMatrixType& X,
               ukfMatrixType& Y) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         (X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          X.cols() == 1) );
  assert(Y.rows() == static_cast<unsigned int>(_signal_dim) &&
         (Y.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          Y.cols() == 1) );
  assert(_signal_data);

  const stdVec_t&       gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t            lambdas1;
  diagmat3_t            lambdas2;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    lambdas1.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas1.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas1.diagonal()[2] = std::max(X(5, i), _lambda_min);
    lambdas2.diagonal()[0] = std::max(X(9, i), _lambda_min);
    lambdas2.diagonal()[1] = std::max(X(10, i), _lambda_min);
    lambdas2.diagonal()[2] = std::max(X(11, i), _lambda_min);

    // Calculate diffusion matrix.
    // HACK: TODO :%s/mat33_t *D/const mat33_t \&D/g
    const mat33_t & D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas1);
    const mat33_t & D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), lambdas2);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j,
        i) =
        std::exp(-b[j] * u.dot(D1 * u) ) * weights_on_tensors_[0] + std::exp(-b[j] * u.dot(D2 * u) )
        * weights_on_tensors_[1];
      }
    }
}

void Full2T::State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1,
                            vec3_t& l1, vec3_t& m2, vec3_t& l2)
{
  // First orientation.
  m1 = rotation_main_dir(x[0], x[1], x[2]);

  // Flip orientation if necessary. (For m1 it should not happen, maybe for
  // m2.)
  if( m1[0] * old_m[0] + m1[1] * old_m[1] + m1[2] * old_m[2] < 0 )
    {
    m1 = -m1;
    }

  // Clamp lambdas.
  l1[0] = std::max(x[3], _lambda_min);
  l1[1] = std::max(x[4], _lambda_min);
  l1[2] = std::max(x[5], _lambda_min);

  // Second orientation.
  m2 = rotation_main_dir(x[6], x[7], x[8]);

  // Flip orientation if necessary.
  if( m2[0] * old_m[0] + m2[1] * old_m[1] + m2[2] * old_m[2] < 0 )
    {
    m2 = -m2;
    }

  // Clamp lambdas.
  l2[0] = std::max(x[9], _lambda_min);
  l2[1] = std::max(x[10], _lambda_min);
  l2[2] = std::max(x[11], _lambda_min);
}
