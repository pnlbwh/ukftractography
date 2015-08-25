#include "filter_Simple2T_FW.h"

// Functions for 2-tensor simple model.
void Simple2T_FW::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
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

    X(3, i) = CheckZero(X(3, i) );
    X(4, i) = CheckZero(X(4, i) );

    X(8, i) = CheckZero(X(8, i) );
    X(9, i) = CheckZero(X(9, i) );

    X(10, i) = CheckZero(X(10, i) );

    if( X(3, i) < 0 || X(4, i) < 0 || X(8, i) < 0 || X(9, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative 1: "
                << X(3, i) << " " << X(4, i) << " " << X(8, i) << " " << X(9, i) << std::endl;
      }
    }
}

void Simple2T_FW::H(const   ukfMatrixType& X,
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
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i) );
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i) );

    lambdas1.diagonal()[0] = CheckZero(X(3, i) );
    lambdas1.diagonal()[1] = CheckZero(X(4, i) );
    lambdas1.diagonal()[2] = lambdas1.diagonal()[1];

    lambdas2.diagonal()[0] = CheckZero(X(8, i) );
    lambdas2.diagonal()[1] = CheckZero(X(9, i) );
    lambdas2.diagonal()[2] = lambdas2.diagonal()[1];

    if( lambdas1.diagonal()[0] < 0 || lambdas1.diagonal()[1] < 0 || lambdas2.diagonal()[0] < 0 ||
        lambdas2.diagonal()[1] < 0 )
      {
      std::cout << "Warning: eigenvalues became negative 2: "
                << lambdas1.diagonal()[0] << " " << lambdas1.diagonal()[1] << " "
                << lambdas2.diagonal()[0] << " " << lambdas2.diagonal()[1] << std::endl;
      }

    // Flip if necessary.
    if( m1[0] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] < 0 )
      {
      m2 = -m2;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(10, i) );

    // Calculate diffusion matrix.
    const mat33_t & D1 = diffusion(m1, lambdas1);
    const mat33_t & D2 = diffusion(m2, lambdas2);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j,
        i) = w
        * (std::exp(-b[j]
                    * u.dot(D1 * u) ) * weights_on_tensors_[0] + std::exp(-b[j] * u.dot(D2 * u) )
           * weights_on_tensors_[1])
        + (1 - w) * std::exp(-b[j] * u.dot(m_D_iso * u) );
      }
    }

}

void Simple2T_FW::State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1,
                                 vec3_t& l1, vec3_t& m2, vec3_t& l2)
{
  // Orientations;
  initNormalized(m1, x[0], x[1], x[2]);
  initNormalized(m2, x[5], x[6], x[7]);

  l1[0] = CheckZero(x[3]);
  l1[1] = CheckZero(x[4]);
  l1[2] = l1[1];

  CheckZero(x[9]);
  l2[0] = CheckZero(x[8]);
  l2[1] = CheckZero(x[9]);
  l2[2] = l2[1];

  if( l1[0] < 0 || l1[1] < 0 || l2[0] < 0 || l2[1] < 0 )
    {
    std::cout << "Warning: eigenvalues became negative 3: " << l1[0] << " " << l1[1] << " " << l2[0] << " "
              << l2[1] << std::endl;
    }

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
}
