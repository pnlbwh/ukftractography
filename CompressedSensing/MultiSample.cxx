#include "MultiSample.h"
#include <cmath>


//function [u] = spiralsample(flg,N)
void spiralsample(unsigned flg, unsigned long N, MatrixType &u)
{
  const double pi(std::atan(static_cast<double>(1.0))*4);
// % SPIRALSAMPLE Spiral sampling of unit hemisphere
// %
// % written by Oleg Michailovich on 12/30/2009
//
// if (nargin==1),
//     N=flg;
//     flg=2;
// end
//
// N=2*N;
  N *= 2;
// z=(1-1/N:-2/N:1/N-1)';
  MatrixType z(N,1);
  double NInv = 1.0/static_cast<double>(N);
  double NIncr = -2.0/static_cast<double>(N);
  double iD = 1.0 - NInv;
  for(unsigned int i = 0; i < N; ++i, iD += NIncr)
    {
    z(i) = iD;
    }
// r=sqrt(1-z.*z);
  MatrixType r(N,1);
  for(unsigned int i = 0; i < N; ++i)
    {
    r(i) = std::sqrt(1.0 - z(i) * z(i));
    }

// switch flg,
  MatrixType Long(N,1);
  switch(flg)
    {
    //     case 1,
    case 1:
      {
      double sqrtN = std::sqrt(static_cast<double>(N));
      // Long=[0; cumsum((3.6/sqrt(N))./r(1:N-1))];
      Long(0) = 0.0;
      for(unsigned i = 1; i < N; ++i)
        {
        Long(i) = Long(i-1) + ((3.6/sqrtN)/r(i-1));
        }
      }
      break;
    case 2:
      {
      // Long=(pi*(3-sqrt(5)))*(0:N-1)';
      double tmp = pi * (3 - std::sqrt(5.0));
      for(unsigned int i = 0; i < N; ++i)
        {
        Long(i) = tmp * static_cast<double>(i);
        }
      }
      break;
    default:
      std::cerr << "Bad Sampling Option" << std::endl;
      throw;
    }
  u.conservativeResize(N,3);
  for(unsigned int i = 0; i < N; ++i)
    {
    // u=[r.*cos(Long) r.*sin(Long) z];
    u(i,0) = r(i) * std::cos(Long(i));
    u(i,1) = r(i) * std::sin(Long(i));
    u(i,2) = z(i);
    // u=u./repmat(sqrt(sum(u.^2,2)),[1 3]);
    u.row(i).normalize();
    }
// u=u(1:N/2,:);
  u = u.block(0,0,N/2,3);
}


extern void MultiSample(unsigned J,
                        std::vector<MatrixType> &v,
                        unsigned long &M,
                        int m0)
{
// function [v,M] = multisample(J,m0)
//
// if (nargin<2),
//     m0=4;
// end
//
// N=(2.^(0:J)*m0+1).^2;
  std::vector<unsigned long> N(J+1);

  for(unsigned int i = 0; i < J+1; ++i)
    {
    unsigned long tmp = std::pow(static_cast<double>(2.0),
                                 static_cast<double>(i)) * m0 + 1;
    N[i] = tmp * tmp;
    }
// v=cell(J+1,1);
  v.resize(J+1);
// for j=1:J+1,
  for(unsigned j = 0; j < J+1; ++j)
    {
//     v{j}=spiralsample(2,N(j));
    spiralsample(2,N[j],v[j]);
// end
    }
// M=sum(N);
  M = 0;
  for(unsigned i = 0; i < J+1; ++i)
    {
    M += N[i];
    }
}

// function [u] = spiralsample(flg,N)
//
// % SPIRALSAMPLE Spiral sampling of unit hemisphere
// %
// % written by Oleg Michailovich on 12/30/2009
//
// if (nargin==1),
//     N=flg;
//     flg=2;
// end
//
// N=2*N;
// z=(1-1/N:-2/N:1/N-1)';
// r=sqrt(1-z.*z);
// switch flg,
//     case 1,
//         long=[0; cumsum((3.6/sqrt(N))./r(1:N-1))];
//     case 2,
//         long=(pi*(3-sqrt(5)))*(0:N-1)';
//     otherwise,
//         error('Invalid sampling option!');
// end
// u=[r.*cos(long) r.*sin(long) z];
// u=u./repmat(sqrt(sum(u.^2,2)),[1 3]);
// u=u(1:N/2,:);

// function [v,M] = multisample(J,m0)
//
// if (nargin<2),
//     m0=4;
// end
//
// N=(2.^(0:J)*m0+1).^2;
// %HACK
// %N(1)=12^2;
// %N(2)=16^2;
// %
// v=cell(J+1,1);
// for j=1:J+1,
//     v{j}=spiralsample(2,N(j));
// end
// M=sum(N);
