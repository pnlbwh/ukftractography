#include "filter_Simple3T.h"

// Functions for 3-tensor simple model.
void Simple3T::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  // Clamp lambdas.
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vectors.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;

    norm_inv = ukfZero; // 1e-16;
    norm_inv += X(5, i) * X(5, i);
    norm_inv += X(6, i) * X(6, i);
    norm_inv += X(7, i) * X(7, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(5, i) *= norm_inv;
    X(6, i) *= norm_inv;
    X(7, i) *= norm_inv;

    norm_inv = ukfZero; // 1e-16;
    norm_inv += X(10, i) * X(10, i);
    norm_inv += X(11, i) * X(11, i);
    norm_inv += X(12, i) * X(12, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(10, i) *= norm_inv;
    X(11, i) *= norm_inv;
    X(12, i) *= norm_inv;

    // Clamp lambdas.
    X(3, i) = std::max(X(3, i), _lambda_min);
    X(4, i) = std::max(X(4, i), _lambda_min);

    X(8, i) = std::max(X(8, i), _lambda_min);
    X(9, i) = std::max(X(9, i), _lambda_min);

    X(13, i) = std::max(X(13, i), _lambda_min);
    X(14, i) = std::max(X(14, i), _lambda_min);
    }
}

void Simple3T::H(const  ukfMatrixType& X,
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

  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  diagmat3_t lambdas3;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i) );
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i) );
    vec3_t m3;
    initNormalized(m3, X(10, i), X(11, i), X(12, i) );

    // Clamp lambdas.
    lambdas1.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas1.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas1.diagonal()[2] = lambdas1.diagonal()[1];

    lambdas2.diagonal()[0] = std::max(X(8, i), _lambda_min);
    lambdas2.diagonal()[1] = std::max(X(9, i), _lambda_min);
    lambdas2.diagonal()[2] = lambdas2.diagonal()[1];

    lambdas3.diagonal()[0] = std::max(X(13, i), _lambda_min);
    lambdas3.diagonal()[1] = std::max(X(14, i), _lambda_min);
    lambdas3.diagonal()[2] = lambdas3.diagonal()[1];

    // flip if necessary
    if( m1[0] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] < 0 )
      {
      m2 = -m2;
      }
    if( m3[0] < 0 )
      {
      m3 = -m3;
      }

    // Calculate diffusion matrix.
    const mat33_t & D1 = diffusion(m1, lambdas1);
    const mat33_t & D2 = diffusion(m2, lambdas2);
    const mat33_t & D3 = diffusion(m3, lambdas3);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) =  std::exp(-b[j] * u.dot(D1 * u) ) * weights_on_tensors_[0]
        + std::exp(-b[j] * u.dot(D2 * u) ) * weights_on_tensors_[1]
        + std::exp(-b[j] * u.dot(D3 * u) ) * weights_on_tensors_[2];
      }
    }
}

void Simple3T::State2Tensor3T(const State& x, const vec3_t& old_m, vec3_t& m1,
                              vec3_t& l1, vec3_t& m2, vec3_t& l2, vec3_t& m3,
                              vec3_t& l3)
{
  // Orientations;
  initNormalized(m1, x[0], x[1], x[2]);
  initNormalized(m2, x[5], x[6], x[7]);
  initNormalized(m3, x[10], x[11], x[12]);

  // Clamp lambdas.
  l1[0] = std::max(x[3], _lambda_min);
  l1[1] = std::max(x[4], _lambda_min);
  l1[2] = l1[1];
  l2[0] = std::max(x[8], _lambda_min);
  l2[1] = std::max(x[9], _lambda_min);
  l2[2] = l2[1];
  l3[0] = std::max(x[13], _lambda_min);
  l3[1] = std::max(x[14], _lambda_min);
  l3[2] = l3[1];

  // Flip orientations if necessary. (For m1 it should not happen, maybe for
  // m2.)
  if( m1[0] * old_m[0] + m1[1] * old_m[1] + m1[2] * old_m[2] < 0 )
    {
    m1 = -m1;
    }
  if( m2[0] * old_m[0] + m2[1] * old_m[1] + m2[2] * old_m[2] < 0 )
    {
    m2 = -m2;
    }
  if( m3[0] * old_m[0] + m3[1] * old_m[1] + m3[2] * old_m[2] < 0 )
    {
    m3 = -m3;
    }
}
