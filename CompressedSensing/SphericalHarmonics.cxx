#include "SphericalHarmonics.h"
#include <cmath>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/geometry/geometry.hpp>
#include <vector>
#include <complex>

MatrixType SphericalHarmonics(const MatrixType &u, unsigned  m)
{
  const double Pi((std::atan(static_cast<double>(1.0))*4));
  const double Pi_2(Pi*2.0);
  const double C(std::sqrt(2.0));
  MatrixType rval(u.rows(),(m+1)*(m+1));

  // i = row index
  for(unsigned i = 0; i < u.rows(); ++i)
    {
    // convert to polar coordinates
    double theta = std::acos(u(i,2));
    double phi = std::atan2(u(i,1),u(i,0));
    if(phi < 0.0)
      {
      phi += Pi_2;
      }
    unsigned int col(0);
    // n = order index
    for(unsigned int n = 0; n < m; n++)
      {
      std::vector<std::complex<double> > curHarm;
      // mm = harmonic index
      for(int mm = -static_cast<int>(n); mm <= static_cast<int>(n); ++mm)
        {
        curHarm.push_back(boost::math::spherical_harmonic(n,mm,theta,phi));
        }
      // the following is my translation of this matlab mumbo-jumbo
      // Y(:,indx)=real([c.*real(y(:,1:k)), y(:,k+1), c.*imag(y(:,k+2:2*k+1))]);
      unsigned int k = 0;
      for(; k < n; ++k,++col)
        {
        rval(i,col) = C * curHarm[k].real();
        }

      rval(i,col) = curHarm[k].real();
      ++col; ++k;

      for(; k < (2 * n) + 1; ++k, ++col)
        {
        rval(i,col) = C * curHarm[k].imag();
        }
      }
    }
  return rval;
}
