/**
 * \file utilities.cc
 * \brief implementation of utilities.h
*/

#include "utilities.h"
#include "math_utilities.h"

double l2fa(double l1, double l2, double l3)
{
  if( l2 == l3 )
    {
    return fabs(l1 - l2) / sqrt(l1 * l1 + 2.0 * l2 * l2);
    }
  else
    {
    return sqrt(0.5 * ( (l1 - l2) * (l1 - l2)
                        + (l2 - l3) * (l2 - l3)
                        + (l3 - l1) * (l3 - l1) )
                / (l1 * l1 + l2 * l2 + l3 * l3) );
    }
}

double s2ga(const vnl_matrix<double>& signal)
{

  int n = signal.rows();

  assert(signal.cols() == 1);

  double mu = 0.0;
  for( int i = 0; i < n; ++i )
    {
    mu += signal(i, 0);
    }
  // average
  mu = mu / n;

  double mu_sq = 0.0;
  double mu_sub = 0.0;
  for( int i = 0; i < n; ++i )
    {
    mu_sq += signal(i, 0) * signal(i, 0);
    mu_sub += (signal(i, 0) - mu) * (signal(i, 0) - mu);
    }

  return sqrt(mu_sub * n) / sqrt( (n - 1) * mu_sq);
}

double curve_radius(const std::vector<vec_t>& fiber)
{
  int length = fiber.size();

  if( length < 3 )
    {
    return 1.0;
    }

  vec_t v1 = fiber[length - 2] - fiber[length - 3];
  vec_t v2 = fiber[length - 1] - fiber[length - 2];

  // Normalize
  v1.normalize();
  v2.normalize();
  double n1 = v1.norm();
  double n2 = v2.norm();
  v1 /= n1;
  v2 /= n2;

  double curv = ( (v2 - v1) / (n2 + n1) ).norm();
  if( isnan(curv) )
    {
    return 0.0;
    }

  return 1.0 / curv;
}
