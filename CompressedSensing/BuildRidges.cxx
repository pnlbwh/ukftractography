#include "CompressedSensing.h"
#include "BuildRidges.h"
#include <cmath>

MatrixType
BuildRidges(unsigned J, double rho, unsigned nu, double p, unsigned flag)
{
  const double pi(std::atan(static_cast<double>(1.0))*4);
  const double pi_4(pi * 4.0);
  // a=2^(1/nu);
  const double a = std::pow(2.0,(1.0/static_cast<double>(nu)));
  // e=1e-9;
  const double localEps(1.0e-9);
  // nmax=a^J*(-log(e)/rho)^(1/p);
  // if rem(ceil(nmax),2)==0,
  //     m=ceil(nmax);
  // else
  //     m=floor(nmax);
  // end
  double nmax = std::pow(a,static_cast<double>(J)) *
    std::pow(-std::log(localEps)/rho,(1.0/p));
  long nmaxInt = static_cast<long>(std::ceil(nmax));
  long m = static_cast<long>(std::floor(nmax));
  if(nmaxInt % 2 == 0)
    {
    m = nmaxInt;
    }
  //
  // n=(0:m)';
  // h=exp(-rho.*(n*a.^(-(0:J) ) ).^p);
  MatrixType h(m+1,3);
  for(unsigned int n = 0; n <= m; ++n)
    {
    const double nD = static_cast<double>(n);
    for(unsigned int j = 0; j < 3; ++j)
      {
      const double jD = static_cast<double>(j);
        h(n,j) = std::exp( -rho * ( nD * std::pow( std::pow(a,jD), p ) ) );
      }
    }
//
// switch flg,
//     case 0,
//         w=[h(:,1) diff(h,[],2)];
//     case 1,
//         w=[h(:,1) sqrt(diff(h.^2,[],2))];
  MatrixType w(m+1,3);
  if(flag == 0)
    {
    for(unsigned int row = 0; row < m+1; ++row)
      {
      w(row,0) = h(row,0);
      w(row,1) = h(row,1) - h(row,0);
      w(row,2) = h(row,2) - h(row,1);
      }
    }
  else
    {
    for(unsigned int row = 0; row < m+1; ++row)
      {
      w(row,0) = h(row,0);
      const double hsq[3] =
        {
          h(row,0) * h(row,0),
          h(row,1) * h(row,1),
          h(row,2) * h(row,2)
        };
      w(row,1) = std::sqrt(hsq[1] - hsq[0]);
      w(row,2) = std::sqrt(hsq[2] - hsq[1]);
      }
    }
  // %w=h;
  // Lmd=ones(m+1,1);
  // Lmd(2:2:m+1)=0;
  // for k=2:2:m,
  //     Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
  // end
  MatrixType Lmd(m+1,1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    Lmd(i)  = (i & 1) == 0 ? 1.0 : 0.0;
    }
  // NOTE: using Matlab indexing, but converting
  // below, since it involves doing math in the index variable.
  for(unsigned k = 2; k <= m; k += 2)
    {
    Lmd(/* k+1 */ k) = -Lmd(/* k-1 */k-2) * (k - 1) / k;
    }
  //
  // psi=w.*Lmd(:,ones(J+1,1));
  //
  MatrixType psi(m+1,J+1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    for(unsigned int j = 0; j < J+1; ++j)
      {
      psi(i,j) = Lmd(i) * w(i,j);
      }
    }
  // n=(0:m)';
  // C=(2*n+1)/(4*pi);
  MatrixType C(m+1,1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    double iD = static_cast<double>(i);
    C(i) = (2.0 * iD + 1.0) / pi_4;
    }
  // psinorm=sqrt(sum(C(:,ones(1,J+1)).*psi.*psi));
  MatrixType psiNorm;
  psiNorm = MatrixType::Zero(1,J+1);
  for(unsigned int row = 0; row < m+1; ++row)
    {
    for(unsigned int col = 0; col < J+1; ++col)
      {
      psiNorm(0,col) += C(row) * psi(row,col) * psi(row,col);
      }
    }
  for(unsigned int col = 0; col < J+1; ++col)
    {
    psiNorm(0,col) = std::sqrt(psiNorm(0,col));
    }
// psi=psi./psinorm(ones(m+1,1),:);
  for(unsigned int row = 0; row < m+1; ++row)
    {
    for(unsigned int col = 0; col < J+1; ++col)
      {
      psi(row,col) = psi(row,col) / psiNorm(0,col);
      }
    }
  return psi;
}

// function [psi,C,m,h,Lmd] = buildridges(J,rho,nu,p,flg,b)
//
// % BUILDRIDGES Legendre coefficients of spherical ridgelets
// %   BUILDRIDGES(J,...) computes the Legendre coefficients of spherical
// % ridgelets which correspond to J+1 resolution levels. The results are
// % stored in an (m+1)-by-(J+1) matrix, where m is defined automatically
// % so as to make the value of the m-th order coefficient at the highest
// % resolution level be bounded by 10^(-9).
// %
// %        function [psi,C,m,h,Lmd] = buildridges(J,rho,nu,p,flg)
// %
// % Inpit:
// %           J - defines the total number of resolution levels as J+1
// %         rho - bandwidth parameter of the ridgelet generating function
// %          nu - number of subresolutions per one octave
// %           p - exponential factor of the ridgelet generating function
// %         flg - defines the type of ridgelet analysis: flg=1 corresponds
// %               to "P-type", while flg=0 corresponds to "M-type".
// % Output:
// %         psi - matrix of Legendre coefficients of spherical ridgelets
// %           C - spherical normalization factor (i.e. C=2*(0:m)/(4*pi))
// %           m - cut-off (truncation) order of the Legendre coefficents
// %           h - Legendre coefficients of the related scaling functions
// %         Lmd - eigenvalues of the Funk-Radon transform
// %
// % see also SPHERICOEFS, POLYLEG, BUILDSENSOR
// %
// % written by Yogesh Rathi, Oleg Michailovich, January 2010.
//
// a=2^(1/nu);
//
// e=1e-9;
// nmax=a^J*(-log(e)/rho)^(1/p);
// if rem(ceil(nmax),2)==0,
//     m=ceil(nmax);
// else
//     m=floor(nmax);
// end
//
// n=(0:m)';
// h=exp(-rho.*(n*a.^(-(0:J))).^p);
//
// switch flg,
//     case 0,
//         w=[h(:,1) diff(h,[],2)];
//     case 1,
//         w=[h(:,1) sqrt(diff(h.^2,[],2))];
//     otherwise
//         error(['The ridgelets type must be equal to either 1 (for "P"',...
//             ' type) or 0 (for "M" type)!']);
// end
// %w=h;
// Lmd=ones(m+1,1);
// Lmd(2:2:m+1)=0;
// for k=2:2:m,
//     Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
// end
//
// psi=w.*Lmd(:,ones(J+1,1));
//
// n=(0:m)';
// C=(2*n+1)/(4*pi);
// psinorm=sqrt(sum(C(:,ones(1,J+1)).*psi.*psi));
// psi=psi./psinorm(ones(m+1,1),:);
