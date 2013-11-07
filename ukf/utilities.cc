/**
 * \file utilities.cc
 * \brief implementation of utilities.h
*/

#include "utilities.h"
#include "math_utilities.h"

ukfPrecisionType l2fa(ukfPrecisionType l1, ukfPrecisionType l2, ukfPrecisionType l3)
{
  if( l2 == l3 )
    {
    return fabs(l1 - l2) / sqrt(l1 * l1 + 2.0 * l2 * l2);
    }
  else
    {
    return sqrt(ukfHalf * ( (l1 - l2) * (l1 - l2)
                        + (l2 - l3) * (l2 - l3)
                        + (l3 - l1) * (l3 - l1) )
                / (l1 * l1 + l2 * l2 + l3 * l3) );
    }
}

ukfPrecisionType s2ga(const ukfMatrixType& signal)
{

  int n = signal.rows();

  assert(signal.cols() == 1);

  ukfPrecisionType mu = ukfZero;
  for( int i = 0; i < n; ++i )
    {
    mu += signal(i, 0);
    }
  // average
  mu = mu / n;

  ukfPrecisionType mu_sq = ukfZero;
  ukfPrecisionType mu_sub = ukfZero;
  for( int i = 0; i < n; ++i )
    {
    mu_sq += signal(i, 0) * signal(i, 0);
    mu_sub += (signal(i, 0) - mu) * (signal(i, 0) - mu);
    }

  return sqrt(mu_sub * n) / sqrt( (n - 1) * mu_sq);
}

ukfPrecisionType curve_radius(const stdVec_t& fiber)
{
  int length = fiber.size();

  if( length < 3 )
    {
    return ukfOne;
    }

  vec3_t v1 = fiber[length - 2] - fiber[length - 3];
  vec3_t v2 = fiber[length - 1] - fiber[length - 2];

  // Normalize
  v1.normalize();
  v2.normalize();
  ukfPrecisionType n1 = v1.norm();
  ukfPrecisionType n2 = v2.norm();

  ukfPrecisionType curv = ( (v2 - v1) / (n2 + n1) ).norm();
  if( std::isnan(curv) )
    {
    return ukfZero;
    }

  return ukfOne / curv;
}
