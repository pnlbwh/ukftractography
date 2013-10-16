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
inline mat_t diag(const double a, const double b, const double c)
{
  return mat_t(vec_t(a,b,c).asDiagonal());
}

/** Assemble rotation matrix given the rotation angles */
inline mat_t rotation(const double theta, const double phi, const double psi)
{
  const double &c_th = cos(theta);
  const double &s_th = sin(theta);
  const double &c_ph = cos(phi);
  const double &s_ph = sin(phi);
  const double &c_ps = cos(psi);
  const double &s_ps = sin(psi);

  const double &q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  const double &q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  const double &q31 = -c_ps * s_th;
  const double &q12 = -c_ps * s_ph - c_th * c_ph * s_ps;
  const double &q22 = c_ph * c_ps - c_th * s_ph * s_ps;
  const double &q32 = s_th * s_ps;
  const double &q13 = c_ph * s_th;
  const double &q23 = s_th * s_ph;
  const double &q33 = c_th;

  mat_t Q;
  Q << q11, q12, q13,
    q21, q22, q23,
    q31, q32, q33;
  return Q;
}

/** Calculate main direction of the diffusion from euler angles */
inline vec_t rotation_main_dir(const double theta, const double phi, const double psi)
{
  const double & c_th = cos(theta);
  const double & s_th = sin(theta);
  const double & c_ph = cos(phi);
  const double & s_ph = sin(phi);
  const double & c_ps = cos(psi);
  const double & s_ps = sin(psi);

  const double & q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  const double & q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  const double & q31 = -c_ps * s_th;
  vec_t rval;
  rval << q11, q21, q31;
  return rval;
}

/** Calculate a diffusion matrix from euler angles */
inline mat_t diffusion_euler(const double theta, const double phi, const double psi,
                             const double l1, const double l2, const double l3)
{
  const mat_t & Q = rotation(theta, phi, psi);
  return Q * diag(l1, l2, l3) * Q.transpose() * 1e-6;
}

/** Make a diffusion tensor matrix from one principal direction, and major and minor EV */
inline mat_t diffusion(const vec_t &m, const double l1, const double l2)
{
  mat_t  R;
  R << m[0], m[1], m[2],
    m[1], m[1] * m[1] / (1 + m[0]) - 1, m[1] * m[2] / (1 + m[0]),
    m[2], m[1] * m[2] / (1 + m[0]), m[2] * m[2] / (1 + m[0]) - 1;

  return R * diag(l1, l2, l2) * R.transpose() * 1e-6;
}

inline void initNormalized(vec_t &m, const double &a, const double &b, const double &c)
{
  m << a, b, c;
  m.normalize();
}

#endif  // LINALG_H_
