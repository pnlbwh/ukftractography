#include "filter_Full1T_FW.h"

// Functions for 1-tensor full model.
void Full1T_FW::F(ukfMatrixType& X) const
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

    X(6, i) = CheckZero(X(6, i) ); // fw

    if( X(3, i) < 0 || X(4, i) < 0 || X(5, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }
    }
}

void Full1T_FW::H(const ukfMatrixType& X,
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
  diagmat3_t            lambdas;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    lambdas.diagonal()[0] = CheckZero(X(3, i) );
    lambdas.diagonal()[1] = CheckZero(X(4, i) );
    lambdas.diagonal()[2] = CheckZero(X(5, i) );

    if( X(3, i) < 0 || X(4, i) < 0 || X(5, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(6, i) );

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) =     (w) * std::exp(-b[j] * u.dot(local_D * u) )
        + (1 - w) * std::exp(-b[j] * u.dot(m_D_iso * u) );
      }
    }
}

void Full1T_FW::State2Tensor1T(const State& x, vec3_t& m, vec3_t& l)
{
  // Orientation.
  m = rotation_main_dir(x[0], x[1], x[2]);

  l[0] = CheckZero(x[3]);
  l[1] = CheckZero(x[4]);
  l[2] = CheckZero(x[5]);
  if( l[0] < 0 || l[1] < 0 || l[2] < 0 )
    {
    std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
    }
}
