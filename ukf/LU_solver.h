/**
 * \file    LU_solver.h
 * \brief   Two algorithms for calculating the LU decomposition of a matrix, and solving a linear system
*/

#ifndef LU_solver_h
#define LU_solver_h

#include <vector>
#define TINY 1.0e-20

/**
 * \namespace LU_Solver
 * \brief Contains functions to solve a linear system AX=B with LU decomposition.
 *
 * This namespace contains two algorithms to solve the linear equation AX=B
 * where size(A)=nxn and size(B)=size(X)=nxm, with n>0,m>0. It was written
 * to efficiently calculate the Gain K in the unscented Kalman Filter of
 * UKF Tractography. As of now, curiously, VNL doesnt provide an efficient solution
 * for matrix RHSs. Be careful with this file, it was only tested with this
 * particular problem.
 * <ol>
 *  <li>Doolittle Algorithm: Slightly faster for the current problem</li>
 *  <li>Crout Algorithm: Numerically more stable because it employs a pivoting
 *      strategy. For this application it made no difference.</li>
 * </ol>
 * \author  Christian Baumgartner (c.f.baumgartner@gmail.com)
*/
namespace LU_Solver
{

/**
 * Calculates the LU decomposition of the square Matrix A : A=LU using the Doolittle algorithm, and
 * stores it in the same Matrix. See http://en.wikipedia.org/wiki/Lu_decomposition#Doolittle_algorithm
 * \param   A nxn matrix to be decomposed and result of LU decomposition (in/out)
 * \param   n Dimension of the matrix.
 * \return  0 if succesful, and -1 if A is singular
*/
int LUdecmpDoolittle(double * const A, const int n)
{
  double * pK = A;

  for( int k = 0; k < n; pK += n, ++k )
    {
    for( int j = k; j < n; ++j )
      {
      const double * pCol = A;
      for( int l = 0; l < k; pCol += n,  ++l )
        {
        *(pK + j) -= *(pK + l) * *(pCol + j);
        }
      }
    if( *(pK + k) == 0.0 )
      {
      return -1;                         // make sure is non-singular
      }
    double * pRow = pK + n;
    for( int i = k + 1; i < n; pRow += n, ++i )
      {
      const double * pCol = A;
      for( int l = 0; l < k; pCol += n, ++l )
        {
        *(pRow + k) -= *(pRow + l) * *(pCol + k);
        }
      *(pRow + k) /= *(pK + k);
      }
    }
  return 0;
}

/**
 * Given the LU decomposition from LUdecmpDoolittle solves the linear system AX=B where
 * X and B are also matrices.
 * \param   LU    the decomposition of A (nRowsxnRows) : A = LU  where the upper triangle contains U,
 *                and the lower L without the diagonal.
 * \param   X     the righthand side matrix B (nRowsxnCols) as input, and X (nRowsxnCols) as output
 * \param   nRows The number of rows in B, X and LU
 * \param   nCols The number of columns in B, and X
 * \return  returns 0 if succesful, and -1 if singular
*/
int LUsolveDoolittle(const double * const LU, double * const X, const int nRows, const int nCols)
{
  for( int i = 0; i < nRows; ++i )
    {
    for( int j = 0; j < nCols; ++j )
      {
      const int iRows = i * nRows;
      const int iCols = i * nCols;
      if( *(LU + iCols + j) == 0.0 )
        {
        return -1;                               // cannot have singular matrix
        }
      for( int k = 0; k < i; ++k )
        {
        *(X + iCols + j) -= *(LU + iRows + k) * *(X + k * nCols + j);
        }
      }
    }
  for( int i = nRows - 1; i >= 0; --i )
    {
    for( int j = nCols - 1; j >= 0; --j )
      {
      const int iRows = i * nRows;
      const int iCols = i * nCols;
      for( int k = nRows - 1; k > i; k-- )
        {
        *(X + iCols + j) -= *(LU + iRows + k) * *(X + k * nCols + j);
        }
      *(X + iCols + j) /= *(LU + iRows + i);
      }
    }
  return 0;
}

/**
 * Calculates the LU decomposition of the square Matrix A : p(A) = LU, where p is a row permutation function.
 * For this version the crout algorithm is used. It was adapted from "Numerical Recipes in C"
 * \param   A     nxn matrix to be decomposed and result of LU decomposition (in/out)
 * \param   n     Dimension of the matrix.
 * \param   indx  A vector stores the original order of the rows, and is used to unscramble a later.
 * \return  0 if succesful, and -1 if A is singular
*/

int LUdecmpCrout(double * const A, const int n, int * const indx)
{
  for( int i = 0; i < n; ++i )
    {
    *(indx + i) = i;
    }

  std::vector<double> vv(n);
  double              big = 0.0;
  for( int i = 0; i < n; ++i )
    {
    for( int j = 0; j < n; ++j )
      {
      const double temp = fabs(*(A + i * n + j) );
      if( temp > big )
        {
        big = temp;
        }
      }
    if( big == 0.0 )
      {
      return -1;
      }
    vv[i] = 1.0 / big;
    }

  int imax = 0;
  for( int j = 0; j < n; ++j )
    {
    for( int i = 0; i < j; ++i )
      {
      double sum = *(A + i * n + j);
      for( int k = 0; k < i; ++k )
        {
        sum -= *(A + i * n + k) * *(A + k * n + j);
        }
      *(A + i * n + j) = sum;
      }
    big = 0.0;
    for( int i = j; i < n; ++i )
      {
      double sum = *(A + i * n + j);
      for( int k = 0; k < j; ++k )
        {
        sum -= *(A + i * n + k) * *(A + k * n + j);
        }
      *(A + i * n + j) = sum;
      const double dum = vv[i] * fabs(sum);
      if( dum >= big )
        {
        big = dum;
        imax = i;
        }
      }
    if( j != imax )
      {
      for( int k = 0; k < n; ++k )
        {
        std::swap((*(A + imax * n + k)),*(A + j * n + k));
        }

      std::swap((*(indx + j)),(*(indx + imax)));

      vv[imax] = vv[j];
      }
    if( *(A + j * n + j) == 0.0 )
      {
      *(A + j * n + j) = TINY;
      }

    if( j != n - 1 )
      {
      const double dum = 1.0 / *(A + j * n + j);
      for( int i = j + 1; i < n; ++i )
        {
        *(A + i * n + j) *= dum;
        }
      }
    }
  return 0;
}

/**
 * Given the LU decomposition from LUdecmpCrout solves the linear system AX=B where
 * X and B are also matrices.
 * \param   LU    the decomposition of A (nRowsxnRows) : A = LU  where the upper triangle contains U,
 *                and the lower L without the diagonal.
 * \param   X     the righthand side matrix B (nRowsxnCols) as input, and X (nRowsxnCols) as output
 * \param   nRows The number of rows in B, X and LU
 * \param   nCols The number of columns in B, and X
 * \param   order The original order of the rows. It's used to undo the permutation of the LU decomposition later. (out)
 * \return  returns 0 if succesful, and -1 if singular
*/
int LUsolveCrout(const double * const LU, const double *const B, double *const X, const int nRows, const int nCols,
                 const int *const order)
{
  for( int i = 0; i < nRows; ++i )
    {
    for( int j = 0; j < nCols; ++j )
      {
      const int iRows = i * nRows;
      const int iCols = i * nCols;
      if( *(LU + iRows + i) == 0.0 )
        {
        return -1;   // i don't like singular matrices
        }
      *(X + iCols + j) = *(B + *(order + i) * nCols + j);
      for( int k = 0; k < i; ++k )
        {
        *(X + iCols + j) -= *(LU + iRows + k) * *(X + k * nCols + j);
        }
      }
    }
  for( int i = nRows - 1; i >= 0; i-- )
    {
    for( int j = nCols - 1; j >= 0; j-- )
      {
      const int iRows = i * nRows;
      const int iCols = i * nCols;
      for( int k = nRows - 1; k > i; k-- )
        {
        *(X + iCols + j) -= *(LU + iRows + k) * *(X + k * nCols + j);
        }
      *(X + iCols + j) /= *(LU + iRows + i);
      }
    }
  return 0;
}

}
#endif
