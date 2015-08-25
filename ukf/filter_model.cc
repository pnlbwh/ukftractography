/**
 * \file filter_model.cc
 * \brief Implements functions defined in filter_model.h
*/

#include "filter_model.h"
#include <iostream>
#include "utilities.h"

ukfPrecisionType FilterModel::CheckZero(const ukfPrecisionType & local_d) const
{
  if( local_d < 0 )
    {
    if( local_d >= -1.0e-4 ) // for small errors just round it to 0
      {
      return ukfZero;
      }
    else   // for errors too big exit with exception
      {
      std::cout << "Error, a variable became negative. Most likely something went wrong in the QP\n";
      exit(1);
      }
    }
  return local_d;
}

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
  const stdVec_t&  gradients = _signal_data->gradients();

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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
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
    //HACK: TODO :%s/mat33_t *D/const mat33_t \&D/g
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas1);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), lambdas2);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j,
        i) =
        std::exp(-b[j] * u.dot(D1 * u) ) * weights_on_tensors_[0] + std::exp(-b[j] * u.dot(D2 * u) ) * weights_on_tensors_[1];
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

// Functions for 3-tensor full model.
void Full3T::F(ukfMatrixType& X) const
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
    X(15, i) = std::max(X(15, i), _lambda_min);
    X(16, i) = std::max(X(16, i), _lambda_min);
    X(17, i) = std::max(X(17, i), _lambda_min);
    }
}

void Full3T::H(const  ukfMatrixType& X,
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  diagmat3_t lambdas3;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    lambdas1.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas1.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas1.diagonal()[2] = std::max(X(5, i), _lambda_min);
    lambdas2.diagonal()[0] = std::max(X(9, i), _lambda_min);
    lambdas2.diagonal()[1] = std::max(X(10, i), _lambda_min);
    lambdas2.diagonal()[2] = std::max(X(11, i), _lambda_min);
    lambdas3.diagonal()[0] = std::max(X(15, i), _lambda_min);
    lambdas3.diagonal()[1] = std::max(X(16, i), _lambda_min);
    lambdas3.diagonal()[2] = std::max(X(17, i), _lambda_min);

    // Calculate diffusion matrix.
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas1);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), lambdas2);
    const mat33_t &D3 = diffusion_euler(X(12, i), X(13, i), X(14, i), lambdas3);
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

void Full3T::State2Tensor3T(const State& x, const vec3_t& old_m, vec3_t& m1,
                          vec3_t& l1, vec3_t& m2, vec3_t& l2, vec3_t& m3,
                          vec3_t& l3)
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

  // Third orientation.
  m3 = rotation_main_dir(x[12], x[13], x[14]);

  // Flip orientation if necessary.
  if( m3[0] * old_m[0] + m3[1] * old_m[1] + m3[2] * old_m[2] < 0 )
    {
    m3 = -m3;
    }

  // Clamp lambdas.
  l3[0] = std::max(x[15], _lambda_min);
  l3[1] = std::max(x[16], _lambda_min);
  l3[2] = std::max(x[17], _lambda_min);
}

// Functions for 1-tensor simple model.
void Simple1T::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vector.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;

    // Clamp lambdas.
    X(3, i) = std::max(X(3, i), _lambda_min);
    X(4, i) = std::max(X(4, i), _lambda_min);
    }
}

void Simple1T::H(const  ukfMatrixType& X,
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize direction.
    vec3_t m;
    initNormalized(m, X(0, i), X(1, i), X(2, i));

    // Clamp lambdas.
    lambdas.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas.diagonal()[2] = lambdas.diagonal()[1];

    // Flip if necessary.
    // Why is that???
    if( m[0] < 0 )
      {
      m = -m;
      }

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion(m, lambdas); // l3 == l2
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) = std::exp(-b[j] * u.dot(local_D * u) ) * weights_on_tensors_[0];
      }
    }
}

void Simple1T::State2Tensor1T(const State& x, vec3_t& m, vec3_t& l)
{
  // Orientation.
  m << x[0], x[1], x[2];
  m.normalize();

  // Clamp lambdas.
  l[0] = std::max(x[3], _lambda_min);
  l[1] = std::max(x[4], _lambda_min);
  l[2] = l[1];
}

// Functions for 2-tensor simple model.
void Simple2T::F(ukfMatrixType& X) const
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

    // Clamp lambdas.
    X(3, i) = std::max(X(3, i), _lambda_min);
    X(4, i) = std::max(X(4, i), _lambda_min);

    X(8, i) = std::max(X(8, i), _lambda_min);
    X(9, i) = std::max(X(9, i), _lambda_min);
    }
}

void Simple2T::H(const  ukfMatrixType& X,
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i));
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i));

    // Clamp lambdas.
    lambdas1.diagonal()[0] = std::max(X(3, i), _lambda_min);
    lambdas1.diagonal()[1] = std::max(X(4, i), _lambda_min);
    lambdas1.diagonal()[2] = lambdas1.diagonal()[1];

    lambdas2.diagonal()[0] = std::max(X(8, i), _lambda_min);
    lambdas2.diagonal()[1] = std::max(X(9, i), _lambda_min);
    lambdas2.diagonal()[2] = lambdas2.diagonal()[1];

    // Flip if necessary.
    if( m1[0] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] < 0 )
      {
      m2 = -m2;
      }

    // Calculate diffusion matrix.
    const mat33_t &D1 = diffusion(m1, lambdas1);
    const mat33_t &D2 = diffusion(m2, lambdas2);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) = std::exp(-b[j] * u.dot(D1 * u) ) * weights_on_tensors_[0]
        + std::exp(-b[j] * u.dot(D2 * u) ) * weights_on_tensors_[1];
      }
    }
}

void Simple2T::State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1,
                            vec3_t& l1, vec3_t& m2, vec3_t& l2)
{
  // Orientations;
  initNormalized(m1,x[0], x[1], x[2]);
  initNormalized(m2,x[5], x[6], x[7]);

  // Clamp lambdas.
  l1[0] = std::max(x[3], _lambda_min);
  l1[1] = std::max(x[4], _lambda_min);
  l1[2] = l1[1];
  l2[0] = std::max(x[8], _lambda_min);
  l2[1] = std::max(x[9], _lambda_min);
  l2[2] = l2[1];

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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();

  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  diagmat3_t lambdas3;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1,X(0, i), X(1, i), X(2, i));
    vec3_t m2;
    initNormalized(m2,X(5, i), X(6, i), X(7, i));
    vec3_t m3;
    initNormalized(m3,X(10, i), X(11, i), X(12, i));

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
    const mat33_t &D1 = diffusion(m1, lambdas1);
    const mat33_t &D2 = diffusion(m2, lambdas2);
    const mat33_t &D3 = diffusion(m3, lambdas3);
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
  initNormalized(m1,x[0], x[1], x[2]);
  initNormalized(m2,x[5], x[6], x[7]);
  initNormalized(m3,x[10], x[11], x[12]);

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

//////// FREE WATER STUFF //////////////////////////////////////////////////////////////
///////  1T SIMPLE MODEL ////

// Functions for 1-tensor simple model.
void Simple1T_FW::F(ukfMatrixType& X) const
{

  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vector.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;

    // Clamp lambdas.
//     X(3, i) = std::max(X(3, i), ukfZero); // 0 because for free water lambdas are constrained to be > 0
//     X(4, i) = std::max(X(4, i), ukfZero); // because imprecission of qp they are sometimes slightly negative on the order
// of e-16
    X(3, i) = CheckZero(X(3, i) );
    X(4, i) = CheckZero(X(4, i) );

    X(5, i) = CheckZero(X(5, i) ); // fw

    // DEBUGGING
    if( X(3, i) < 0 || X(4, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }
    }
}

void Simple1T_FW::H(const ukfMatrixType& X,
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b     = _signal_data->GetBValues();
  diagmat3_t lambdas;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize direction.
    vec3_t m;
    initNormalized(m,X(0, i), X(1, i), X(2, i));

    // Clamp lambdas.
    lambdas.diagonal()[0]= CheckZero(X(3, i) );
    lambdas.diagonal()[1]= CheckZero(X(4, i) );
    lambdas.diagonal()[2]= lambdas.diagonal()[1];
    if( lambdas.diagonal()[0] < 0 ||  lambdas.diagonal()[1]  < 0 )
      {
      std::cout << "Warning: eigenvalues became negative l1= " << lambdas.diagonal()[0] << "l2= " << lambdas.diagonal()[1] << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero( X(5, i) );
    // FOR DEBUGGING :
    if( w < 0 - 1.0e-5 )
      {
      std::cout << "Negative weight! w= " << w << "\n";
      std::cout << X << "\ni: " << i;
      exit(1);
      }
    if( w > 1 + 1.0e-5 )
      {
      std::cout << "Weight > 1 => negative free water! w= " << w << "\n";
      std::cout << X << "\ni: " << i;
      exit(1);
      }

    // Flip if necessary.
    // Why is that???
    if( m[0] < 0 )
      {
      m = -m;
      }

    // Calculate diffusion matrix.
    const mat33_t &local_D = diffusion(m, lambdas); // l3 == l2
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) =  (w) * std::exp(-b[j] * u.dot(local_D * u) )
        + (1 - w) * std::exp(-b[j] * u.dot(m_D_iso * u) );
      }
    }
}

void Simple1T_FW::State2Tensor1T(const State& x, vec3_t& m, vec3_t& l)
{
  // Orientation.
  initNormalized(m,x[0], x[1], x[2]);

  l[1] = std::max(x[4], ukfZero );

  l[0] = CheckZero(x[3]);
  l[1] = CheckZero(x[4]);
  if( l[0] < 0 || l[1] < 0 )
    {
    std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
    }
  l[2] = l[1];
}

///////  1T FULL MODEL ///
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas;
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
    const mat33_t &local_D = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas);
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

////////// 2T SIMPLE MODEL ///
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

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1,X(0, i), X(1, i), X(2, i) );
    vec3_t m2;
    initNormalized(m2,X(5, i), X(6, i), X(7, i) );

    lambdas1.diagonal()[0] = CheckZero(X(3, i) );
    lambdas1.diagonal()[1] = CheckZero(X(4, i) );
    lambdas1.diagonal()[2] = lambdas1.diagonal()[1];

    lambdas2.diagonal()[0] = CheckZero(X(8, i) );
    lambdas2.diagonal()[1] = CheckZero(X(9, i) );
    lambdas2.diagonal()[2] = lambdas2.diagonal()[1];

    if( lambdas1.diagonal()[0] < 0 || lambdas1.diagonal()[1] < 0 || lambdas2.diagonal()[0] < 0 || lambdas2.diagonal()[1] < 0 )
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
    const mat33_t &D1 = diffusion(m1, lambdas1);
    const mat33_t &D2 = diffusion(m2, lambdas2);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j,
        i) = w
        * (std::exp(-b[j]
               * u.dot(D1 * u) ) * weights_on_tensors_[0] + std::exp(-b[j] * u.dot(D2 * u) ) * weights_on_tensors_[1])
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

  // Flip orientations if necessary. (For m1 it should not happen, maybe for m2.)
  if( m1[0] * old_m[0] + m1[1] * old_m[1] + m1[2] * old_m[2] < 0 )
    {
    m1 = -m1;
    }
  if( m2[0] * old_m[0] + m2[1] * old_m[1] + m2[2] * old_m[2] < 0 )
    {
    m2 = -m2;
    }
}

////////// 2T FULL MODEL ///
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
          X.cols() == 1));
  assert(Y.rows() == static_cast<unsigned int>(_signal_dim) &&
         (Y.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          Y.cols() == 1) );
  assert(_signal_data);

  const stdVec_t&   gradients = _signal_data->gradients();
  const ukfVectorType & b       = _signal_data->GetBValues();
  diagmat3_t lambdas1;
  diagmat3_t lambdas2;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {

    lambdas1.diagonal()[0]= CheckZero(X(3, i) );
    lambdas1.diagonal()[1]= CheckZero(X(4, i) );
    lambdas1.diagonal()[2]= CheckZero(X(5, i) );
    lambdas2.diagonal()[0]= CheckZero(X(9, i) );
    lambdas2.diagonal()[1]= CheckZero(X(10, i) );
    lambdas2.diagonal()[2]= CheckZero(X(11, i) );

    if( X(3, i) < 0 || X(4, i) < 0 ||  X(5, i) < 0 || X(9, i) < 0 || X(10, i) < 0 ||  X(11, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(12, i) );

    // Calculate diffusion matrix.
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), lambdas1);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), lambdas2);

    // Reconstruct signal.
    ukfMatrixType valMatrix(_signal_dim,3);
    for( int j = 0; j < _signal_dim; ++j )
    {
       const vec3_t& u = gradients[j];
       valMatrix(j,0)=-b[j] * u.dot(D1 * u);
       valMatrix(j,1)=-b[j] * u.dot(D2 * u);
       valMatrix(j,2)=-b[j] * u.dot(m_D_iso * u);
    }
    const ukfMatrixType expMatrix = valMatrix.array().exp();
    for( int j = 0; j < _signal_dim; ++j )
      {
      const ukfPrecisionType part1a = expMatrix(j,0) * weights_on_tensors_[0];
      const ukfPrecisionType part1b = expMatrix(j,1) * weights_on_tensors_[1];
      const ukfPrecisionType part1 = part1a + part1b;
      const ukfPrecisionType part2 = expMatrix(j,2);

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

void createProtocol(const ukfVectorType& _b_values,
                    ukfVectorType& _gradientStrength, ukfVectorType& _pulseSeparation)
{
  static int count=0;
  static ukfVectorType gradientStrength, pulseSeparation;
  std::vector<double> Bunique, tmpG;
  ukfPrecisionType Bmax = 0;
  ukfPrecisionType tmp, Gmax, GAMMA;

  _gradientStrength.resize(_b_values.size());
  _pulseSeparation.resize(_b_values.size());

  // set maximum G = 40 mT/m
  Gmax = 0.04;
  GAMMA = 267598700;
  if(count ==1)
  {
     // gradient strength and pulseSeparation are once
     // same values are returned
    _gradientStrength = gradientStrength;
    _pulseSeparation = pulseSeparation;
  }
  else{
    for(int i = 0; i < _b_values.size(); ++i )
    {
      int unique = 1;
      for(size_t j = 0; j < Bunique.size(); ++j )
      {
        if (_b_values[i] == Bunique[j])
        {
          unique = 0;
          break;
        }
      }
      if(unique == 1)
      {
        Bunique.push_back(_b_values[i]);
      }
      if(Bmax < _b_values[i])
      {
        Bmax = _b_values[i];
      }
    }

    tmp = cbrt(3*Bmax*1000000/(2*GAMMA*GAMMA*Gmax*Gmax));

    for(int i = 0; i < _b_values.size(); ++i )
    {
      _pulseSeparation[i] = tmp;
    }

    for(size_t i = 0; i < Bunique.size(); ++i )
    {
      tmpG.push_back(std::sqrt(Bunique[i]/Bmax) * Gmax);
    }

    for(size_t i = 0; i < Bunique.size(); ++i )
    {
      for(int j=0; j < _b_values.size(); j++)
      {
        if(_b_values[j] == Bunique[i])
        {
          _gradientStrength[j] = tmpG[i];
        }
      }
    }
    gradientStrength = _gradientStrength;
    pulseSeparation = _pulseSeparation;
    count = 1;
  }
}

// Functions for 1-fiber NODDI model.
void NODDI1F::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vector.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;
    }
}

void NODDI1F::H(const  ukfMatrixType& X, ukfMatrixType& Y) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         (X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          X.cols() == 1) );
  assert(Y.rows() == static_cast<unsigned int>(_signal_dim) &&
         (Y.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          Y.cols() == 1) );
  assert(_signal_data);

  ukfPrecisionType dPar, dIso;
  dPar = 0.0000000017;
  dIso = 0.000000003;
  const stdVec_t&   gradients = _signal_data->gradients();
  ukfVectorType gradientStrength, pulseSeparation ;
  createProtocol(_signal_data->GetBValues(),gradientStrength,pulseSeparation);
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    ukfVectorType Eec, Eic, Eiso;
    // Normalize direction.
    vec3_t m;
    initNormalized(m, X(0, i), X(1, i), X(2, i));

    // Clamp lambdas.
    const ukfPrecisionType Vic =X(3, i);
    const ukfPrecisionType kappa = X(4, i);
    const ukfPrecisionType Viso =X(5, i);
    ExtraCelluarModel(dPar, Vic, kappa, gradientStrength,
                          pulseSeparation, gradients, m, Eec);
    IntraCelluarModel(dPar, kappa, gradientStrength, pulseSeparation,
                          gradients, m, Eic);
    IsoModel(dIso, gradientStrength, pulseSeparation, Eiso);
    assert(Eic.size() == _signal_dim);
    for( int j = 0; j < _signal_dim; ++j )
      {
      Y(j, i) = Viso*Eiso[j] + (1- Viso)*(Vic*Eic[j] + (1 - Vic)*Eec[j]);
      }
    }
}

// Functions for 2-fiber noddi model.
void NODDI2F::F(ukfMatrixType& X) const
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
    }
}

void NODDI2F::H(const  ukfMatrixType& X,
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

  ukfPrecisionType dPar, dIso;
  dPar = 0.0000000017;
  dIso = 0.000000003;
  const stdVec_t&   gradients = _signal_data->gradients();
  ukfVectorType gradientStrength, pulseSeparation ;
  createProtocol(_signal_data->GetBValues(),gradientStrength,pulseSeparation);
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i));
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i));

    // Flip if necessary.
    if( m1[0] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] < 0 )
      {
      m2 = -m2;
      }

    ukfVectorType Eec1, Eic1, Eec2, Eic2, Eiso;
    const ukfPrecisionType Vic1 = X(3, i);
    const ukfPrecisionType kappa1 = X(4, i);
    const ukfPrecisionType Vic2 = X(8, i);
    const ukfPrecisionType kappa2 = X(9, i);
    const ukfPrecisionType Viso = X(10, i);

    ExtraCelluarModel(dPar, Vic1, kappa1, gradientStrength,
                          pulseSeparation, gradients, m1, Eec1);
    IntraCelluarModel(dPar, kappa1, gradientStrength, pulseSeparation,
                          gradients, m1, Eic1 );

    ExtraCelluarModel(dPar, Vic2, kappa2, gradientStrength,
                          pulseSeparation, gradients, m2, Eec2);
    IntraCelluarModel(dPar, kappa2, gradientStrength, pulseSeparation,
                          gradients, m2, Eic2 );
    IsoModel(dIso, gradientStrength, pulseSeparation, Eiso);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      Y(j, i) = Viso*Eiso[j]+(1-Viso)*((Vic1*Eic1[j] + (1 - Vic1)*Eec1[j]) * weights_on_tensors_[0] +
                         (Vic2*Eic2[j] + (1 - Vic2)*Eec2[j]) * weights_on_tensors_[1]);
      }
    }
}
