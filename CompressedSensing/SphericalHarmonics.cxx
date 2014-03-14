#include "SphericalHarmonics.h"
#include <cmath>

#ifdef BOOST_SPHERICALHARMONIC
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
      std::vector<unsigned long> ind;

      std::vector<std::complex<double> > curHarm;
      // mm = harmonic index
      for(int mm = -static_cast<int>(n); mm <= static_cast<int>(n); ++mm)
        {
        curHarm.push_back(boost::math::spherical_harmonic(n,mm,theta,phi));
        }
      // the following is my translation of this matlab mumbo-jumbo
      // Y(:,indx)=real([c.*real(y(:,1:k)), y(:,k+1), c.*imag(y(:,k+2:2*k+1))]);
      const unsigned start(n*n);
      const unsigned limit((n +1 ) * (n + 1));
      for(col = start; col < limit; ++col)
        {
        const unsigned k(col - start);
        if(k < n)
          {
          rval(i,col) = C * curHarm[k].real();
          }
        else if(k == n)
          {
          rval(i,col) = curHarm[k].real();
          }
        else
          {
          rval(i,col) = C * curHarm[k].imag();
          }
        }
      }
    }
  return rval;
}
#endif
#ifdef FINSLER_SPHERICALHARMONIC
#include "SphHarm.h"

MatrixType SphericalHarmonics(const MatrixType &u, unsigned  m)
{
  MatrixType rval;
  unsigned int rows = u.rows();
  double *x = new double[rows];
  double *y = new double[rows];
  double *z = new double[rows];
  double *theta = new double[rows];
  double *phi = new double[rows];
  for(unsigned int i = 0; i < rows; ++i)
    {
    x[i] = u(i,0);
    y[i] = u(i,1);
    z[i] = u(i,2);
    }
  shmaths::computeSphericalCoordsFromCartesian(x,y,z,theta,phi,rows);
  shmaths::computeSHMatrix(rows,theta,phi,m,rval);
  delete [] x; delete [] y; delete [] z;
  delete [] theta; delete [] phi;
  return rval;
}
#endif
#define MATLAB_SPHERICALHARMONIC
#ifdef MATLAB_SPHERICALHARMONIC

#include "SphericalHarmonicsConst.h"

MatrixType SphericalHarmonics(const MatrixType &u, unsigned  m)
{
  if(m != 22)
    {
    std::cerr << "Have to rebuild SphericalHarmonicsConst.h from Matlab" << std::endl;
    exit(1);
    }
  unsigned long cols = (m+1)*(m+1);
  MatrixType rval(u.rows(),cols);
  for(unsigned int i = 0; i < u.rows(); ++i)
    {
    for(unsigned int j = 0; j < cols; ++j)
      {
      rval(i,j) = sphericalHarmonicsConstant[i][j];
      }
    }
  return rval;
}
#endif
