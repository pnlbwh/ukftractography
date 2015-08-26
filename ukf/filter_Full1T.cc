#include "filter_Full1T.h"

// Functions for 1-tensor full model.
void Full1T::F(ukfMatrixType& X) const
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
    }
}

void Full1T::H(const  ukfMatrixType& X,
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

  const ukfVectorType& b        = _signal_data->GetBValues();
  const stdVec_t&      gradients = _signal_data->gradients();

  diagmat3_t lambdas;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    lambdas.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas.diagonal()[2] = std::max(X(5, i), _lambda_min);

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) = std::exp(-b[j] * u.dot(local_D * u) ) * weights_on_tensors_[0];
      }
    }
}

void Full1T::State2Tensor1T(const State& x, vec3_t& m, vec3_t& l)
{
  // Orientation.
  m = rotation_main_dir(x[0], x[1], x[2]);

  // Clamp lambdas.
  l[0] = std::max(x[3], _lambda_min);
  l[1] = std::max(x[4], _lambda_min);
  l[2] = std::max(x[5], _lambda_min);
}
