/**
 * \file filter_model.cc
 * \brief Implements functions defined in filter_model.h
*/

#include "filter_model.h"
#include <iostream>

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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    const ukfPrecisionType l1 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l2 = std::max(X(4, i), _lambda_min);
    const ukfPrecisionType l3 = std::max(X(5, i), _lambda_min);

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion_euler(X(0, i), X(1, i), X(2, i), l1, l2, l3);
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    const ukfPrecisionType l11 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l12 = std::max(X(4, i), _lambda_min);
    const ukfPrecisionType l13 = std::max(X(5, i), _lambda_min);
    const ukfPrecisionType l21 = std::max(X(9, i), _lambda_min);
    const ukfPrecisionType l22 = std::max(X(10, i), _lambda_min);
    const ukfPrecisionType l23 = std::max(X(11, i), _lambda_min);

    // Calculate diffusion matrix.
    //HACK: TODO :%s/mat33_t *D/const mat33_t \&D/g
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), l11, l12, l13);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), l21, l22, l23);
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Clamp lambdas.
    const ukfPrecisionType l11 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l12 = std::max(X(4, i), _lambda_min);
    const ukfPrecisionType l13 = std::max(X(5, i), _lambda_min);
    const ukfPrecisionType l21 = std::max(X(9, i), _lambda_min);
    const ukfPrecisionType l22 = std::max(X(10, i), _lambda_min);
    const ukfPrecisionType l23 = std::max(X(11, i), _lambda_min);
    const ukfPrecisionType l31 = std::max(X(15, i), _lambda_min);
    const ukfPrecisionType l32 = std::max(X(16, i), _lambda_min);
    const ukfPrecisionType l33 = std::max(X(17, i), _lambda_min);

    // Calculate diffusion matrix.
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), l11, l12, l13);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), l21, l22, l23);
    const mat33_t &D3 = diffusion_euler(X(12, i), X(13, i), X(14, i), l31, l32, l33);
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize direction.
    vec3_t m;
    initNormalized(m, X(0, i), X(1, i), X(2, i));

    // Clamp lambdas.
    const ukfPrecisionType l1 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l2 = std::max(X(4, i), _lambda_min);

    // Flip if necessary.
    // Why is that???
    if( m[0] < 0 )
      {
      m = -m;
      }

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion(m, l1, l2); // l3 == l2
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

  ukfPrecisionType n = m.norm();
  m /= n;

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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i));
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i));

    // Clamp lambdas.
    const ukfPrecisionType l11 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l12 = std::max(X(4, i), _lambda_min);

    const ukfPrecisionType l21 = std::max(X(8, i), _lambda_min);
    const ukfPrecisionType l22 = std::max(X(9, i), _lambda_min);

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
    const mat33_t &D1 = diffusion(m1, l11, l12);
    const mat33_t &D2 = diffusion(m2, l21, l22);
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
    const ukfPrecisionType l11 = std::max(X(3, i), _lambda_min);
    const ukfPrecisionType l12 = std::max(X(4, i), _lambda_min);

    const ukfPrecisionType l21 = std::max(X(8, i), _lambda_min);
    const ukfPrecisionType l22 = std::max(X(9, i), _lambda_min);

    const ukfPrecisionType l31 = std::max(X(13, i), _lambda_min);
    const ukfPrecisionType l32 = std::max(X(14, i), _lambda_min);

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
    const mat33_t &D1 = diffusion(m1, l11, l12);
    const mat33_t &D2 = diffusion(m2, l21, l22);
    const mat33_t &D3 = diffusion(m3, l31, l32);
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
      std::cout << "Warning: eigenvalues became negative" << std::endl;
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize direction.
    vec3_t m;
    initNormalized(m,X(0, i), X(1, i), X(2, i));

    // Clamp lambdas.
    //     const ukfPrecisionType l1 = std::max(X(3, i), ukfZero);
    //     const ukfPrecisionType l2 = std::max(X(4, i), ukfZero);

    const ukfPrecisionType l1 = CheckZero(X(3, i) );
    const ukfPrecisionType l2 = CheckZero(X(4, i) );
    if( l1 < 0 || l2 < 0 )
      {
      std::cout << "Warning: eigenvalues became negative l1= " << l1 << "l2= " << l2  << std::endl;
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
    const mat33_t &local_D = diffusion(m, l1, l2); // l3 == l2
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
    std::cout << "Warning: eigenvalues became negative" << std::endl;
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
      std::cout << "Warning: eigenvalues became negative" << std::endl;
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {

    const ukfPrecisionType l1 = CheckZero(X(3, i) );
    const ukfPrecisionType l2 = CheckZero(X(4, i) );
    const ukfPrecisionType l3 = CheckZero(X(5, i) );

    if( X(3, i) < 0 || X(4, i) < 0 || X(5, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative" << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(6, i) );

    // Calculate diffusion matrix.
    const mat33_t &local_D = diffusion_euler(X(0, i), X(1, i), X(2, i), l1, l2, l3);
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
    std::cout << "Warning: eigenvalues became negative" << std::endl;
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1,X(0, i), X(1, i), X(2, i) );
    vec3_t m2;
    initNormalized(m2,X(5, i), X(6, i), X(7, i) );

    const ukfPrecisionType l11 = CheckZero(X(3, i) );
    const ukfPrecisionType l12 = CheckZero(X(4, i) );

    const ukfPrecisionType l21 = CheckZero(X(8, i) );
    const ukfPrecisionType l22 = CheckZero(X(9, i) );

    if( l11 < 0 || l12 < 0 || l21 < 0 || l22 < 0 )
      {
      std::cout << "Warning: eigenvalues became negative 2: " << l11 << " " << l12 << " " << l21 << " " << l22
                << std::endl;
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
    const mat33_t &D1 = diffusion(m1, l11, l12);
    const mat33_t &D2 = diffusion(m2, l21, l22);
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
      std::cout << "Warning: eigenvalues became negative" << std::endl;
      }
    }
}

void Full2T_FW::H(const   ukfMatrixType& X,
                  ukfMatrixType& Y) const
{
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
  for( unsigned int i = 0; i < X.cols(); ++i )
    {

    const ukfPrecisionType l11 = CheckZero(X(3, i) );
    const ukfPrecisionType l12 = CheckZero(X(4, i) );
    const ukfPrecisionType l13 = CheckZero(X(5, i) );
    const ukfPrecisionType l21 = CheckZero(X(9, i) );
    const ukfPrecisionType l22 = CheckZero(X(10, i) );
    const ukfPrecisionType l23 = CheckZero(X(11, i) );

    if( X(3, i) < 0 || X(4, i) < 0 ||  X(5, i) < 0 || X(9, i) < 0 || X(10, i) < 0 ||  X(11, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative" << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero(X(12, i) );

    // Calculate diffusion matrix.
    const mat33_t &D1 = diffusion_euler(X(0, i), X(1, i), X(2, i), l11, l12, l13);
    const mat33_t &D2 = diffusion_euler(X(6, i), X(7, i), X(8, i), l21, l22, l23);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) = (w)
        * (std::exp(-b[j] * u.dot(D1 * u) )
           * weights_on_tensors_[0] + std::exp(-b[j] * u.dot(D2 * u) )
           * weights_on_tensors_[1])
        + (1 - w) * std::exp(-b[j] * u.dot(m_D_iso * u) );
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
    std::cout << "Warning: eigenvalues became negative" << std::endl;
    }
}
