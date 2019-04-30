/**
 * \file linalg.h
 *
 * \brief Contains commonly used definitions and functions for linear algebra in 3 dimensions
*/

#ifndef LINALG_H_
#define LINALG_H_

#include <cmath>
#include "ukf_types.h"

/** Assemble rotation matrix given the rotation angles */
inline mat33_t rotation(const ukfPrecisionType theta, const ukfPrecisionType phi, const ukfPrecisionType psi)
{
  const ukfPrecisionType &c_th = std::cos(theta);
  const ukfPrecisionType &s_th = std::sin(theta);
  const ukfPrecisionType &c_ph = std::cos(phi);
  const ukfPrecisionType &s_ph = std::sin(phi);
  const ukfPrecisionType &c_ps = std::cos(psi);
  const ukfPrecisionType &s_ps = std::sin(psi);

  const ukfPrecisionType &q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  const ukfPrecisionType &q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  const ukfPrecisionType &q31 = -c_ps * s_th;
  const ukfPrecisionType &q12 = -c_ps * s_ph - c_th * c_ph * s_ps;
  const ukfPrecisionType &q22 = c_ph * c_ps - c_th * s_ph * s_ps;
  const ukfPrecisionType &q32 = s_th * s_ps;
  const ukfPrecisionType &q13 = c_ph * s_th;
  const ukfPrecisionType &q23 = s_th * s_ph;
  const ukfPrecisionType &q33 = c_th;

  mat33_t Q;
  Q << q11, q12, q13,
    q21, q22, q23,
    q31, q32, q33;
  return Q;
}

/** Calculate main direction of the diffusion from euler angles */
inline vec3_t rotation_main_dir(const ukfPrecisionType theta, const ukfPrecisionType phi, const ukfPrecisionType psi)
{
  const ukfPrecisionType & c_th = std::cos(theta);
  const ukfPrecisionType & s_th = std::sin(theta);
  const ukfPrecisionType & c_ph = std::cos(phi);
  const ukfPrecisionType & s_ph = std::sin(phi);
  const ukfPrecisionType & c_ps = std::cos(psi);
  const ukfPrecisionType & s_ps = std::sin(psi);

  const ukfPrecisionType & q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  const ukfPrecisionType & q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  const ukfPrecisionType & q31 = -c_ps * s_th;
  vec3_t rval;
  rval << q11, q21, q31;
  return rval;
}

/** Calculate a diffusion matrix from euler angles */
inline mat33_t diffusion_euler(const ukfPrecisionType theta, const ukfPrecisionType phi, const ukfPrecisionType psi,
                             const diagmat3_t & lambdas)
{
  const mat33_t & Q = rotation(theta, phi, psi);
  return Q * lambdas * Q.transpose() * GLOBAL_TENSOR_UNPACK_VALUE; //NOTE: Scale factor to offset Tensor::UnpackTensor scaling
}

/** Make a diffusion tensor matrix from one principal direction, and major and minor EV */
inline mat33_t diffusion(const vec3_t &m, const diagmat3_t & lambdas)
{
  mat33_t  R;
  R << m[0], m[1], m[2],
    m[1], m[1] * m[1] / (1 + m[0]) - 1, m[1] * m[2] / (1 + m[0]),
    m[2], m[1] * m[2] / (1 + m[0]), m[2] * m[2] / (1 + m[0]) - 1;

  return R * lambdas * R.transpose() * GLOBAL_TENSOR_UNPACK_VALUE;
}

inline mat33_t diffusion_l2eql3(const vec3_t &eigenVec1, const ukfPrecisionType L1, const ukfPrecisionType L2)
{
  // $ D = ( ( \lambda_1 - \lambda_2) * e_1 *e_1^{T} + \lambda_2 * I $
  const mat33_t & D = ((L1-L2)* eigenVec1 * eigenVec1.transpose() + L2 * mat33_t::Identity())*GLOBAL_TENSOR_UNPACK_VALUE;
  return D;
}

inline void initNormalized(vec3_t &m, const ukfPrecisionType &a, const ukfPrecisionType &b, const ukfPrecisionType &c)
{
  m << a, b, c;
  m.normalize();
}

#endif  // LINALG_H_
