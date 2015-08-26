#include "filter_Full2T_FW.h"

void Full2T_FW::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  // Clamp lambdas.
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    X(3, i) = CheckZero(X(3, i) );
    X(4, i) = CheckZero(X(4, i) );
    X(5, i) = CheckZero(X(5, i) );
    X(9, i) = CheckZero(X(9, i) );
    X(10, i) = CheckZero(X(10, i) );
    X(11, i) = CheckZero(X(11, i) );

    X(12, i) = CheckZero(X(12, i) );

    if( X(3, i) < 0 || X(4, i) < 0 ||  X(5, i) < 0 || X(9, i) < 0 || X(10, i) < 0 ||  X(11, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }
    }
}

extern unsigned int countH;
void Full2T_FW::H(const   ukfMatrixType& X,
                  ukfMatrixType& Y) const
{
  countH++;
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

    lambdas1.diagonal()[0] = CheckZero(X(3, i) );
    lambdas1.diagonal()[1] = CheckZero(X(4, i) );
    lambdas1.diagonal()[2] = CheckZero(X(5, i) );
    lambdas2.diagonal()[0] = CheckZero(X(9, i) );
    lambdas2.diagonal()[1] = CheckZero(X(10, i) );
    lambdas2.diagonal()[2] = CheckZero(X(11, i) );

    if( X(3, i) < 0 || X(4, i) < 0 ||  X(5, i) < 0 || X(9, i) < 0 || X(10, i) < 0 ||  X(11, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(12, i) );

    // Calculate diffusion matrix.
    const mat33_t & D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas1);
    const mat33_t & D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), lambdas2);

    // Reconstruct signal.
    ukfMatrixType valMatrix(_signal_dim, 3);
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      valMatrix(j, 0) = -b[j] * u.dot(D1 * u);
      valMatrix(j, 1) = -b[j] * u.dot(D2 * u);
      valMatrix(j, 2) = -b[j] * u.dot(m_D_iso * u);
      }
    const ukfMatrixType expMatrix = valMatrix.array().exp();
    for( int j = 0; j < _signal_dim; ++j )
      {
      const ukfPrecisionType part1a = expMatrix(j, 0) * weights_on_tensors_[0];
      const ukfPrecisionType part1b = expMatrix(j, 1) * weights_on_tensors_[1];
      const ukfPrecisionType part1 = part1a + part1b;
      const ukfPrecisionType part2 = expMatrix(j, 2);

      Y(j, i) = (w) * part1 + ( ukfOne - w) * part2;
      }
    }
}

void Full2T_FW::State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1,
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

  l1[0] = CheckZero(x[3]);
  l1[1] = CheckZero(x[4]);
  l1[2] = CheckZero(x[5]);

  // Second orientation.
  m2 = rotation_main_dir(x[6], x[7], x[8]);

  // Flip orientation if necessary.
  if( m2[0] * old_m[0] + m2[1] * old_m[1] + m2[2] * old_m[2] < 0 )
    {
    m2 = -m2;
    }

  l2[0] = CheckZero(x[9]);
  l2[1] = CheckZero(x[10]);
  l2[2] = CheckZero(x[11]);

  if( l1[0] < 0 || l1[1] < 0 ||  l1[2] < 0 || l2[0] < 0 || l2[1] < 0 ||  l2[2] < 0 )
    {
    std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
    }
}
