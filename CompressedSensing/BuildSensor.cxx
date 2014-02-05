#include "BuildSensor.h"
#include <cmath>
MatrixType polyleg(MatrixType x,unsigned n)
{
// x=x(:); // this seems to do nothing
// N=length(x);
  unsigned int N = x.rows();
//
  MatrixType absX = x;
  for(unsigned int i = 0; i < N; ++i)
    {
    if(absX(i) < 0.0)
      {
      absX(i) *= -1.0;
      }
    }
// if (max(abs(x))>1),
  if(absX.maxCoeff() > 1.0)
    {
//     x(abs(x)>1) = sign(x(abs(x)>1));
//     %warning('The argument is required to be in [-1,1]!');
// end
    for(unsigned int i = 0; i < N; ++i)
      {
      if(absX(i) > 1.0)
        {
        x(i) = x(i) > 0.0 ? 1.0 : -1.0;
        }
      }
    }
//
// P=ones(N,n+1);
  MatrixType P;
  P = MatrixType::Ones(N,n+1);
// switch n,
  switch(n)
    {
    // case 0,
    case 0:
      // NOTE: this doesn't do anything
      // P(:,1)=ones(N,1);
      P.col(0) = MatrixType::Ones(N,1);
      break;
    // case 1,
    case 1:
      // P(:,2)=x;
      P.col(1) = x;
      break;
    // otherwise,
    default:
      // P(:,2)=x;
      P.col(1) = x;
      // for k=3:n+1,
      for(unsigned k = 3; k <= n+1; ++k)
        {
        // c1=(2*k-3)/(k-1);
        const double kD(k);
        double c1 = (2.0 * (kD - 3)) / (kD - 1.0);
        // c2=(k-2)/(k-1);
        double c2 = (kD - 2.0) / (kD - 1.0);
        // P(:,k)=c1.*x.*P(:,k-1)-c2*P(:,k-2);
        P.col(k - 1) = c1 * x * P.col(k-2) - c2 * P.col(k-3);
        // for(unsigned int i = 0; i < N; ++i)
        //   {
        //   P(i,k-1) = c1 * x(i) * P(i,k-2) - c2 * P(i,k -3);
        //   }
        // end
        }
    }
  return P;
}
// function [P] = polyleg(x,n)
// % POLYLEG - the Legendre polynomials
// %   P = POLYLEG(x,n) is a length(x)-by-(n+1) matrix the columns of which
// %   contain the Legendre polynomials of orders 0,1,...,n evaluated at x.
// %   Note that x must obey |x|<=1.
// %
// % written by Oleg Michailovich, 08/22/2007.
//
// x=x(:);
// N=length(x);
//
// if (max(abs(x))>1),
//     x(abs(x)>1) = sign(x(abs(x)>1));
//     %warning('The argument is required to be in [-1,1]!');
// end
//
// P=ones(N,n+1);
// switch n,
//     case 0,
//         P(:,1)=ones(N,1);
//     case 1,
//         P(:,2)=x;
//     otherwise,
//         P(:,2)=x;
//         for k=3:n+1,
//             c1=(2*k-3)/(k-1);
//             c2=(k-2)/(k-1);
//             P(:,k)=c1.*x.*P(:,k-1)-c2*P(:,k-2);
//         end
// end

MatrixType BuildSensor(MatrixType &g, std::vector<MatrixType> &v,MatrixType &psi)
{
  const double pi(std::atan(static_cast<double>(1.0))*4);
  const double pi_4(pi * 4.0);
  // K=size(g,1);
  unsigned K = g.rows();
  // J=length(v)-1;
  unsigned J = v.size() - 1;
  // m=size(psi,1)-1;
  unsigned m = psi.rows() - 1;

// C=(2*(0:m)'+1)/(4*pi);
  MatrixType C(m+1,1);
  for(unsigned i = 0; i < m+1; ++i)
    {
    C(i) = (2.0 * static_cast<double>(i) + 1.0) / pi_4;
    }
  // XCn=C(:,ones(1,J+1)).*psi;
  MatrixType XCn(m+1,J+1);
  for(unsigned int i = 0; i < m+1; ++i)
    {
    for(unsigned int j = 0; j < J+1; ++j)
      {
      XCn(i,j) = C(i) * psi(i,j);
      }
    }
  // N=zeros(J+1,1);
  // for k=0:J,
  //     N(k+1)=size(v{k+1},1);
  // end
  MatrixType N;
  unsigned nSum = 0;
  N = MatrixType::Zero(J+1,1);
  for(unsigned int k = 0; k < J+1; ++k)
    {
    nSum += v[k].size();
    N(k) = v[k].size();

    }
  // **** As it happens, there's only ever one return value requested
  // switch nargout,
  //     case 1,
  //         A=zeros(K,sum(N));
  //         for j=0:J,
  //             vv=v{j+1};
  //             for k=1:N(j+1),
  //                 P=polyleg(g*vv(k,:)',m);
  //                 A(:,k+sum(N(1:j)))=P*XCn(:,j+1);
  //             end
  //         end
  MatrixType A;
  A = MatrixType::Zero(K,nSum);
  for(unsigned int j = 0; j < J+1; ++j)
    {
    MatrixType &vv = v[j];
    for(unsigned int k = 0; k < N(j); ++k)
      {
      MatrixType tmp1 = g * vv.row(k).transpose();
      MatrixType P = polyleg(tmp1,m);
      A.col(nSum) = P * XCn.col(j);
      }
    }
  return A;
  //     case 2,
  //         Lmd=ones(m+1,1);
  //         Lmd(2:2:m+1)=0;
  //         t=(0:m).*(1:m+1);t=t';
  //         for k=2:2:m,
  //             Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
  //         end
  //         XCL=(t(:,ones(1,J+1)).*XCn.*Lmd(:,ones(1,J+1)));
  //
  //         A=zeros(K,sum(N));
  //         Q=zeros(K,sum(N));
  //         for j=0:J,
  //             vv=v{j+1};
  //             for k=1:N(j+1),
  //                 P=polyleg(g*vv(k,:)',m);
  //                 A(:,k+sum(N(1:j)))=P*XCn(:,j+1);
  //                 Q(:,k+sum(N(1:j)))=P*XCL(:,j+1);
  //             end
  //         end
  // end
}
// function [A,Q] = buildsensor(g,v,psi)
//
// % BUILDSENSOR Sensing matrix for ridgelets-based compressed sensing
// %  BUILDSENSOR(v,g,psi) computes the values of spherical ridgelets defined
// % by the Legendre coefficients psi and orientations v. The values are com-
// % puted at the spherical points defined by an N-by-3 matrix g.
// %
// %                         [A,Q] = buildsensor(g,v,psi)
// % Input:
// %               v - (J+1)-by-1 cell array of spherical coordinates where
// %                   u{j} (j=1,...,J+1) defines the rotation set of sphe-
// %                   rical ridgelets at resolution j.
// %               g - sampling set over which the evaluation is performed
// %             psi - m-by-(J+1) matrix of the Legendre coefficients of the
// %                   spherical ridgelets
// % Output:
// %               A - compressed sensing matrix
// %               Q - Funk-Radon transforms of the ridgelets in A
// %
// % See also BUILDRIDGES, SPHERICOEFS, POLYLEG, SPHERAPPROX
// %
// % writen by Yogesh Rathi, Oleg Michailovich, February 2010.
// %rot3D=@(d)(1/(d'*[0; 0; 1]+1))*(([0; 0; 1]+d)*([0; 0; 1]+d)')-eye(3);
// K=size(g,1);
// J=length(v)-1;
// m=size(psi,1)-1;
//
// C=(2*(0:m)'+1)/(4*pi);
// XCn=C(:,ones(1,J+1)).*psi;
//
// N=zeros(J+1,1);
// for k=0:J,
//     N(k+1)=size(v{k+1},1);
// end
//
// switch nargout,
//     case 1,
//         A=zeros(K,sum(N));
//         for j=0:J,
//             vv=v{j+1};
//             for k=1:N(j+1),
//                 P=polyleg(g*vv(k,:)',m);
//                 A(:,k+sum(N(1:j)))=P*XCn(:,j+1);
//             end
//         end
//     case 2,
//         Lmd=ones(m+1,1);
//         Lmd(2:2:m+1)=0;
//         t=(0:m).*(1:m+1);t=t';
//         for k=2:2:m,
//             Lmd(k+1)=-Lmd(k-1)*(k-1)/k;
//         end
//         XCL=(t(:,ones(1,J+1)).*XCn.*Lmd(:,ones(1,J+1)));
//
//         A=zeros(K,sum(N));
//         Q=zeros(K,sum(N));
//         for j=0:J,
//             vv=v{j+1};
//             for k=1:N(j+1),
//                 P=polyleg(g*vv(k,:)',m);
//                 A(:,k+sum(N(1:j)))=P*XCn(:,j+1);
//                 Q(:,k+sum(N(1:j)))=P*XCL(:,j+1);
//             end
//         end
// end
//
