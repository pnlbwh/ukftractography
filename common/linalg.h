/**
 * \file linalg.h
 *
 * \brief Contains commonly used definitions and functions for linear algebra in 3 dimensions
*/

#ifndef LINALG_H_
#define LINALG_H_

#if defined(WIN32)
#define hypot _hypot
#endif

#include <cmath>
#include "ukf_types.h"

/** Make a diagonal matrix */
inline mat33_t diag(const ukfPrecisionType a, const ukfPrecisionType b, const ukfPrecisionType c)
{
  return mat33_t(vec3_t(a,b,c).asDiagonal());
}

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
                             const ukfPrecisionType l1, const ukfPrecisionType l2, const ukfPrecisionType l3)
{
  const mat33_t & Q = rotation(theta, phi, psi);
  return Q * diag(l1, l2, l3) * Q.transpose() * GLOBAL_TENSOR_UNPACK_VALUE; //NOTE: Scale factor to offset Tensor::UnpackTensor scaling
}

/** Make a diffusion tensor matrix from one principal direction, and major and minor EV */
inline mat33_t diffusion(const vec3_t &m, const ukfPrecisionType l1, const ukfPrecisionType l2)
{
  mat33_t  R;
  R << m[0], m[1], m[2],
    m[1], m[1] * m[1] / (1 + m[0]) - 1, m[1] * m[2] / (1 + m[0]),
    m[2], m[1] * m[2] / (1 + m[0]), m[2] * m[2] / (1 + m[0]) - 1;

  return R * diag(l1, l2, l2) * R.transpose() * GLOBAL_TENSOR_UNPACK_VALUE;
}

inline void initNormalized(vec3_t &m, const ukfPrecisionType &a, const ukfPrecisionType &b, const ukfPrecisionType &c)
{
  m << a, b, c;
  m.normalize();
}

#endif  // LINALG_H_
