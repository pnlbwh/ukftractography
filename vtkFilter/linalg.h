#ifndef LINALG_H_
#define LINALG_H_

#include <cmath>

// Defines simple 3D vector and matrices.
typedef struct {
  double _[3];
} vec_t;

inline vec_t make_vec(double x, double y, double z)
{
  vec_t m = { {x, y, z} };
  return m;
}

inline vec_t make_vec(double th, double phi)
{
  double cosphi = cos(phi);
  double x = cosphi * cos(th), y = cosphi * sin(th), z = sin(phi);
  return make_vec(x, y, z);
}

inline vec_t make_vec_s(double *sph)
{
  double th = sph[0], phi = sph[1];
  return make_vec(th, phi);
}

inline vec_t make_vec_c(double *m)
{
  return make_vec(m[0], m[1], m[2]);
}

inline void vec2sph(vec_t m, double *sph)
{
  double x = m._[0], y = m._[1], z = m._[2];
  sph[0] = atan2(y, x);
  sph[1] = atan2(z, hypot(x, y));
}

inline void vec2mem(vec_t v, double *m)
{
  for (int i = 0; i < 3; i++) {
    m[i] = v._[i];
  }
}

// Negation.
inline vec_t operator-(vec_t a)
{
  return make_vec(-a._[0], -a._[1], -a._[2]);
}

// Addition.
inline vec_t operator+(vec_t a, vec_t b)
{
  return make_vec(a._[0] + b._[0], a._[1] + b._[1], a._[2] + b._[2]);
}

inline vec_t operator+(vec_t a, double b)
{
  return make_vec(a._[0] + b, a._[1] + b, a._[2] + b);
}

inline vec_t operator+(double a, vec_t b)
{
  return b + a;
}

inline void operator+=(vec_t &a, double b)
{
  a._[0] += b;
  a._[1] += b;
  a._[2] += b;
}

// Subtraction.
inline vec_t operator-(vec_t a, vec_t b)
{
  return make_vec(a._[0] - b._[0], a._[1] - b._[1], a._[2] - b._[2]);
}

inline vec_t operator-(vec_t a, double b)
{
  return make_vec(a._[0] - b, a._[1] - b, a._[2] - b);
}

inline vec_t operator-(double a, vec_t b)
{
  return make_vec(a - b._[0], a - b._[1], a - b._[2]);
}

inline void operator-=(vec_t &a, double b)
{
  a._[0] -= b;
  a._[1] -= b;
  a._[2] -= b;
}

// Multiplication.
inline vec_t operator*(vec_t a, double b)
{
  return make_vec(a._[0] * b, a._[1] * b, a._[2] * b);
}

inline vec_t operator*(double a, vec_t b)
{
  return b * a;
}

inline void operator*=(vec_t &a, double b)
{
  a._[0] *= b;
  a._[1] *= b;
  a._[2] *= b;
}

// Division.
inline vec_t operator/(vec_t a, double b)
{
  return make_vec(a._[0] / b, a._[1] / b, a._[2] / b);
}

inline vec_t operator/(double a, vec_t b)
{
  return make_vec(a / b._[0], a / b._[1], a / b._[2]);
}

inline void operator/=(vec_t &a, double b)
{
  a._[0] /= b;
  a._[1] /= b;
  a._[2] /= b;
}

inline double norm(vec_t a)
{
  return sqrt(a._[0] * a._[0] + a._[1] * a._[1] + a._[2] * a._[2]);
}

inline double norm2(vec_t a)
{
  return a._[0] * a._[0] + a._[1] * a._[1] + a._[2] * a._[2];
}

inline double dot(vec_t a, vec_t b)
{
  return a._[0] * b._[0] + a._[1] * b._[1] + a._[2] * b._[2];
}

typedef struct {
  double _[9];
} mat_t;

inline mat_t make_mat(double a, double b, double c,
                      double d, double e, double f,
                      double g, double h, double i)
{
  mat_t m = { {a, b, c, d, e, f, g, h, i} };
  return m;
}

inline mat_t diag(double a, double b, double c)
{
  mat_t m = { {a, 0, 0, 0, b, 0, 0, 0, c} };
  return m;
}

inline mat_t t(mat_t m)
{
  return make_mat(m._[0], m._[3], m._[6],
                  m._[1], m._[4], m._[7],
                  m._[2], m._[5], m._[8]);
}

// Matrix-scalar operations.
inline mat_t operator*(mat_t a, double b)
{
  return make_mat(a._[0] * b,  a._[1] * b, a._[2] * b,
                  a._[3] * b,  a._[4] * b, a._[5] * b,
                  a._[6] * b,  a._[7] * b, a._[8] * b);
}

inline mat_t operator*(double a, mat_t b)
{
  return b * a;
}

// Matrix-vector operations.
inline vec_t operator*(mat_t a, vec_t b)
{
  return make_vec(a._[0] * b._[0] + a._[1] * b._[1] + a._[2] * b._[2],
                  a._[3] * b._[0] + a._[4] * b._[1] + a._[5] * b._[2],
                  a._[6] * b._[0] + a._[7] * b._[1] + a._[8] * b._[2]);
}

// Matrix-matrix operations.
inline mat_t operator*(mat_t a, mat_t b)
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

// Determinant.
inline double det(mat_t M)
{
  return M._[0] * (M._[4] * M._[8] - M._[5] * M._[7]) -
         M._[1] * (M._[3] * M._[8] - M._[5] * M._[6]) +
         M._[2] * (M._[3] * M._[7] - M._[4] * M._[6]);
}

// Conjugate transpose.
inline mat_t ct(mat_t M)
{
  return make_mat( M._[0], -M._[3],  M._[6],
                   -M._[1],  M._[4], -M._[7],
                   M._[2], -M._[5],  M._[8]);
}

// Inversion.
inline mat_t inv(mat_t M)
{
  return (1 / det(M)) * ct(M);
}

inline mat_t diffusion(vec_t m, double l1, double l2)
{
  double x = m._[0], y = m._[1], z = m._[2];
  mat_t R = make_mat(x, y, z,
                     y, y * y / (1 + x) - 1, y * z / (1 + x),
                     z, y * z / (1 + x), z * z / (1 + x) - 1);
  return R * diag(l1, l2, l2) * t(R) * 1e-6;
}

// Calculates a rotation matrix from euler angles.
inline mat_t rotation(double theta, double phi, double psi)
{
  double c_th = cos(theta);
  double s_th = sin(theta);
  double c_ph = cos(phi);
  double s_ph = sin(phi);
  double c_ps = cos(psi);
  double s_ps = sin(psi);

  double q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  double q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  double q31 = -c_ps * s_th;
  double q12 = -c_ps * s_ph - c_th * c_ph * s_ps;
  double q22 = c_ph * c_ps - c_th * s_ph * s_ps;
  double q32 = s_th * s_ps;
  double q13 = c_ph * s_th;
  double q23 = s_th * s_ph;
  double q33 = c_th;

  mat_t Q = make_mat(q11, q12, q13,
                     q21, q22, q23,
                     q31, q32, q33);
  return Q;
}

// Only gives the main direction of the diffusion from euler angles.
inline vec_t rotation_main_dir(double theta, double phi, double psi)
{
  double c_th = cos(theta);
  double s_th = sin(theta);
  double c_ph = cos(phi);
  double s_ph = sin(phi);
  double c_ps = cos(psi);
  double s_ps = sin(psi);

  double q11 = c_th * c_ph * c_ps - s_ph * s_ps;
  double q21 = c_th * c_ps * s_ph + c_ph * s_ps;
  double q31 = -c_ps * s_th;

  return make_vec(q11, q21, q31);
}

// Calculates a diffusion matrix from euler angles.
inline mat_t diffusion_euler(double theta, double phi, double psi,
                             double l1, double l2, double l3)
{
  mat_t Q = rotation(theta, phi, psi);
  return Q * diag(l1, l2, l3) * t(Q) * 1e-6;
}

#endif  // LINALG_H_
