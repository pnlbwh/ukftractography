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
#include "Eigen/Dense"
/**
 * \struct vec_t
 * \brief Defines simple 3D vector and matrices. 3
*/
// typedef struct
//  {
//  double _[3];
//  } vec_t;
typedef Eigen::Vector3d  vec_t;

#if 0
/** Make a 3D vector */
inline vec_t make_vec(const double x, const double y, const double z)
{
  const vec_t m = { {x, y, z} };
  return m;
}

/** Make a 3D vector from two angles */
inline vec_t make_vec(const double th, const double phi)
{
  const double cosphi = cos(phi);
  const double x = cosphi * cos(th);
  const double y = cosphi * sin(th);
  const double z = sin(phi);
  return make_vec(x, y, z);
}

/** Make a 3D vector from a double array of two with angles */
inline vec_t make_vec_s(const double *const sph)
{
  const double th = sph[0];
  const double phi = sph[1];
  return make_vec(th, phi);
}

/** Make a 3D vector from a double array of three */
inline vec_t make_vec_c(const double *const m)
{
  return make_vec(m[0], m[1], m[2]);
}

/** Make a double array of angles from a 3D vector */
inline void vec2sph(const vec_t &m, double *sph)
{
  const double &x = m._[0];
  const double &y = m._[1];
  const double &z = m._[2];

  sph[0] = atan2(y, x);
  sph[1] = atan2(z, hypot(x, y) );
}

/** Make a double array of three from a 3D vector */
inline void vec2mem(const vec_t &v, double * const m)
{
  for( int i = 0; i < 3; ++i )
    {
    m[i] = v._[i];
    }
}

// Operators

/** Negation */
inline vec_t operator-(const vec_t & a)
{
  return make_vec(-a._[0], -a._[1], -a._[2]);
}

/** Addition two vectors*/
inline vec_t operator+(const vec_t &a, const vec_t &b)
{
  return make_vec(a._[0] + b._[0], a._[1] + b._[1], a._[2] + b._[2]);
}

/** Addition, add scalar to each vector element */
inline vec_t operator+(const vec_t &a, const double b)
{
  return make_vec(a._[0] + b, a._[1] + b, a._[2] + b);
}

/** Addition, add scalar to each vector element */
inline vec_t operator+(const double a, const vec_t &b)
{
  return b + a;
}

/** Agumented assignment of scalar to vector */
inline void operator+=(vec_t & a, const double b)
{
  a._[0] += b;
  a._[1] += b;
  a._[2] += b;
}

/** Substraction of two vectors */
inline vec_t operator-(const vec_t &a, const vec_t &b)
{
  return make_vec(a._[0] - b._[0], a._[1] - b._[1], a._[2] - b._[2]);
}

/** Substract a scalar from each vector element */
inline vec_t operator-(const vec_t &a, const double b)
{
  return make_vec(a._[0] - b, a._[1] - b, a._[2] - b);
}

/** Substract each vector element from a scalar */
inline vec_t operator-(const double a, const vec_t &b)
{
  return make_vec(a - b._[0], a - b._[1], a - b._[2]);
}

/** Substraction compound assignment */
inline void operator-=(vec_t & a, const double b)
{
  a._[0] -= b;
  a._[1] -= b;
  a._[2] -= b;
}

/** Scaling of a vector by a scalar*/
inline vec_t operator*(const vec_t &a, const double b)
{
  return make_vec(a._[0] * b, a._[1] * b, a._[2] * b);
}

/** Scaling of a vector by a scalar */
inline vec_t operator*(const double a, const vec_t &b)
{
  return b * a;
}

/** Compound assignment scaling */
inline void operator*=(vec_t & a, const double b)
{
  a._[0] *= b;
  a._[1] *= b;
  a._[2] *= b;
}

/** Scaling of a vector by division */
inline vec_t operator/(const vec_t & a, const double b)
{
  return make_vec(a._[0] / b, a._[1] / b, a._[2] / b);
}

/** Scaling of a vector by division */
inline vec_t operator/(const double a, const vec_t &b)
{
  return make_vec(a / b._[0], a / b._[1], a / b._[2]);
}

/** Compound assignment scaling by division */
inline void operator/=(vec_t & a, const double b)
{
  a._[0] /= b;
  a._[1] /= b;
  a._[2] /= b;
}

/** Calculate the euler-norm of a vector */
inline double norm(const vec_t &a)
{
  return sqrt(a._[0] * a._[0] + a._[1] * a._[1] + a._[2] * a._[2]);
}

/** Calculate the 2-norm of a vector */
inline double norm2(const vec_t &a)
{
  return a._[0] * a._[0] + a._[1] * a._[1] + a._[2] * a._[2];
}

/** Calculate the dot product of two vectors */
inline double dot(const vec_t &a, const vec_t &b)
{
  return a._[0] * b._[0] + a._[1] * b._[1] + a._[2] * b._[2];
}
#endif

/**
 * \struct mat_t
 * \brief Defines a 3-by-3 Matrix
*/
// typedef struct
//   {
//   double _[9];
//   } mat_t;
typedef Eigen::Matrix<double,3,3> mat_t;
#if 0
/** Make a matrix */
inline mat_t make_mat(const double a, const double b, const double c,
                      const double d, const double e, const double f,
                      const double g, const double h, const double i)
{
  const mat_t m = { {a, b, c, d, e, f, g, h, i} };
  return m;
}

/** Make a diagonal matrix */
inline mat_t diag(const double a, const double b, const double c)
{
  const mat_t m = { {a, 0, 0, 0, b, 0, 0, 0, c} };
  return m;
}

/** Transpose a matrix */
inline mat_t t(const mat_t & m)
{
  return make_mat(m._[0], m._[3], m._[6],
                  m._[1], m._[4], m._[7],
                  m._[2], m._[5], m._[8]);
}

/** Scale all elements of a matrix */
inline mat_t operator*(const mat_t & a, const double b)
{
  return make_mat(a._[0] * b,  a._[1] * b, a._[2] * b,
                  a._[3] * b,  a._[4] * b, a._[5] * b,
                  a._[6] * b,  a._[7] * b, a._[8] * b);
}

/** Scale all elements of a matrix */
inline mat_t operator*(const double a, const mat_t & b)
{
  return b * a;
}

/** Multiply a matrix (a) with vector (b) */
inline vec_t operator*(const mat_t & a, const vec_t & b)
{
  return make_vec(a._[0] * b._[0] + a._[1] * b._[1] + a._[2] * b._[2],
                  a._[3] * b._[0] + a._[4] * b._[1] + a._[5] * b._[2],
                  a._[6] * b._[0] + a._[7] * b._[1] + a._[8] * b._[2]);
}

/** Mutliply two 3-by-3 Matrices */
inline mat_t operator*(const mat_t & a, const mat_t & b)
{
  return make_mat(a._[0] * b._[0] + a._[1] * b._[3] + a._[2] * b._[6],
                  a._[0] * b._[1] + a._[1] * b._[4] + a._[2] * b._[7],
                  a._[0] * b._[2] + a._[1] * b._[5] + a._[2] * b._[8],

                  a._[3] * b._[0] + a._[4] * b._[3] + a._[5] * b._[6],
                  a._[3] * b._[1] + a._[4] * b._[4] + a._[5] * b._[7],
                  a._[3] * b._[2] + a._[4] * b._[5] + a._[5] * b._[8],

                  a._[6] * b._[0] + a._[7] * b._[3] + a._[8] * b._[6],
                  a._[6] * b._[1] + a._[7] * b._[4] + a._[8] * b._[7],
                  a._[6] * b._[2] + a._[7] * b._[5] + a._[8] * b._[8]);
}

/** Calculate the determinant of a 3-by-3 matrix */
inline double det(const mat_t & M)
{
  return M._[0] * (M._[4] * M._[8] - M._[5] * M._[7])
         - M._[1] * (M._[3] * M._[8] - M._[5] * M._[6])
         + M._[2] * (M._[3] * M._[7] - M._[4] * M._[6]);
}

/** Conjugate transpose of a 3-by-3 matrix */
inline mat_t ct(const mat_t &M)
{
  return make_mat( M._[0], -M._[3],  M._[6],
                   -M._[1],  M._[4], -M._[7],
                   M._[2], -M._[5],  M._[8]);
}

/** Inversion of a 3-by-3 matrix */
inline mat_t inv(const mat_t &M)
{
  return (1 / det(M) ) * ct(M);
}

/** Make a diffusion tensor matrix from one principal direction, and major and minor EV */
inline mat_t diffusion(const vec_t &m, const double l1, const double l2)
{
  const double &x = m._[0];
  const double &y = m._[1];
  const double &z = m._[2];
  mat_t  R = make_mat(x, y, z,
                      y, y * y / (1 + x) - 1, y * z / (1 + x),
                      z, y * z / (1 + x), z * z / (1 + x) - 1);

  return R * diag(l1, l2, l2) * t(R) * 1e-6;
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

  const mat_t & Q = make_mat(q11, q12, q13,
                     q21, q22, q23,
                     q31, q32, q33);
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

  return make_vec(q11, q21, q31);
}

/** Calculate a diffusion matrix from euler angles */
inline mat_t diffusion_euler(const double theta, const double phi, const double psi,
                             const double l1, const double l2, const double l3)
{
  const mat_t & Q = rotation(theta, phi, psi);
  return Q * diag(l1, l2, l3) * t(Q) * 1e-6;
}
#endif

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
