/**
 * \file QuadProg++_Eigen.cc
 * \brief implementation of QuadProg++_Eigen.h
 *
 * This file was adapted from QuadProg++ an open project available on
 * sourceforge.net (see http://sourceforge.net/projects/quadprog/). The major change
 * is that the file now works entirely with Eigen, and is not dependant on the
 * helper classes Vector and Matrix in Array.hh. Furthermore the equality constraints
 * have been removed. The ce0 and CE variables passed are simply dummy variables.
 * If you need equality constraints change it back in the code. See the bottom of
 * the file for additional comments by the original authors.
 *
 * \author Christian Baumgartner (baumgach@ee.ethz.ch) adapted from code by Luca Di Gaspero
*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include "QuadProg++_Eigen.h"

#include "Eigen/Dense"

// #define TRACE_SOLVER
namespace QuadProgPP
{
// Utility functions for updating some data needed by the solution method
inline void compute_d(ukfVectorType& d, const ukfMatrixType& J, const ukfVectorType& np)
{
    d = J.adjoint() * np;
}

inline void update_z(ukfVectorType& z, const ukfMatrixType& J, const ukfVectorType& d, int iq)
{
    z = J.rightCols(z.size() - iq) * d.tail(d.size() - iq);
}

inline void update_r(const ukfMatrixType& R, ukfVectorType& r, const ukfVectorType& d, int iq)
{
    /* setting of r = R^-1 d */
    r.head(iq) = R.topLeftCorner(iq, iq).triangularView<Eigen::Upper>().solve(d.head(iq));
}

inline ukfPrecisionType distance(ukfPrecisionType a, ukfPrecisionType b)
{
    const ukfPrecisionType a1 = ::std::fabs(a);
    const ukfPrecisionType b1 = ::std::fabs(b);
    if (a1 > b1)
    {
        ukfPrecisionType t = (b1 / a1);
        return a1 * ::std::sqrt(1.0 + t * t);
    }
    else if (b1 > a1)
    {
        ukfPrecisionType t = (a1 / b1);
        return b1 * ::std::sqrt(1.0 + t * t);
    }
    return a1 * ::std::sqrt(2.0);
}

inline bool add_constraint(ukfMatrixType& R, ukfMatrixType& J, ukfVectorType& d, int& iq, ukfPrecisionType& R_norm)
{
    using std::abs;
    int n = J.rows();
#ifdef TRACE_SOLVER
    std::cerr << "Add constraint " << iq << '/';
#endif
    int j, k;
    double cc, ss, h, t1, t2, xny;

    /* we have to find the Givens rotation which will reduce the element
          d(j) to zero.
          if it is already zero we don't have to do anything, except of
          decreasing j */
    for (j = n - 1; j >= iq + 1; j--)
    {
        /* The Givens rotation is done with the matrix (cc cs, cs -cc).
                 If cc is one, then element (j) of d is zero compared with element
                 (j - 1). Hence we don't have to do anything.
                 If cc is zero, then we just have to switch column (j) and column (j - 1)
                 of J. Since we only switch columns in J, we have to be careful how we
                 update d depending on the sign of gs.
                 Otherwise we have to apply the Givens rotation to these columns.
                 The i - 1 element of d has to be updated to h. */
        cc = d(j - 1);
        ss = d(j);
        h = distance(cc, ss);
        if (h == 0.0)
            continue;
        d(j) = 0.0;
        ss = ss / h;
        cc = cc / h;
        if (cc < 0.0)
        {
            cc = -cc;
            ss = -ss;
            d(j - 1) = -h;
        }
        else
            d(j - 1) = h;
        xny = ss / (1.0 + cc);
        for (k = 0; k < n; k++)
        {
            t1 = J(k, j - 1);
            t2 = J(k, j);
            J(k, j - 1) = t1 * cc + t2 * ss;
            J(k, j) = xny * (t1 + J(k, j - 1)) - t2;
        }
    }
    /* update the number of constraints added*/
    iq++;
    /* To update R we have to put the iq components of the d vector
      into column iq - 1 of R
      */
    R.col(iq - 1).head(iq) = d.head(iq);
#ifdef TRACE_SOLVER
    std::cerr << iq << std::endl;
#endif

    if (abs(d(iq - 1)) <= std::numeric_limits<double>::epsilon() * R_norm)
        // problem degenerate
        return false;
    R_norm = std::max<double>(R_norm, abs(d(iq - 1)));
    return true;
}

inline void delete_constraint(ukfMatrixType& R, ukfMatrixType& J, Eigen::VectorXi& A, ukfVectorType& u,
                       int p, int& iq, int l)
{
    int n = R.rows();
#ifdef TRACE_SOLVER
    std::cerr << "Delete constraint " << l << ' ' << iq;
#endif
    int i, j, k;
    int qq = -1; // just to prevent warnings from smart compilers
    double cc, ss, h, xny, t1, t2;

    /* Find the index qq for active constraint l to be removed */
    for (i = p; i < iq; i++)
        if (A(i) == l)
        {
            qq = i;
            break;
        }

    /* remove the constraint from the active set and the duals */
    for (i = qq; i < iq - 1; i++)
    {
        A(i) = A(i + 1);
        u(i) = u(i + 1);
        R.col(i) = R.col(i + 1);
    }

    A(iq - 1) = A(iq);
    u(iq - 1) = u(iq);
    A(iq) = 0;
    u(iq) = 0.0;
    for (j = 0; j < iq; j++)
        R(j, iq - 1) = 0.0;
    /* constraint has been fully removed */
    iq--;
#ifdef TRACE_SOLVER
    std::cerr << '/' << iq << std::endl;
#endif

    if (iq == 0)
        return;

    for (j = qq; j < iq; j++)
    {
        cc = R(j, j);
        ss = R(j + 1, j);
        h = distance(cc, ss);
        if (h == 0.0)
            continue;
        cc = cc / h;
        ss = ss / h;
        R(j + 1, j) = 0.0;
        if (cc < 0.0)
        {
            R(j, j) = -h;
            cc = -cc;
            ss = -ss;
        }
        else
            R(j, j) = h;

        xny = ss / (1.0 + cc);
        for (k = j + 1; k < iq; k++)
        {
            t1 = R(j, k);
            t2 = R(j + 1, k);
            R(j, k) = t1 * cc + t2 * ss;
            R(j + 1, k) = xny * (t1 + R(j, k)) - t2;
        }
        for (k = 0; k < n; k++)
        {
            t1 = J(k, j);
            t2 = J(k, j + 1);
            J(k, j) = t1 * cc + t2 * ss;
            J(k, j + 1) = xny * (J(k, j) + t1) - t2;
        }
    }
}

// The Solving function, implementing the Goldfarb-Idnani method

ukfPrecisionType solve_quadprog(ukfMatrixType& G, ukfVectorType& g0,
                      const ukfMatrixType& CE, const ukfVectorType& ce0,
                      const ukfMatrixType& CI, const ukfVectorType& ci0,
                      ukfVectorType& x)
{
  std::ostringstream msg;

  // Ensure that the dimensions of the matrices and vectors can be
  // safely converted from unsigned int into to int without overflow.
  const unsigned mx = std::numeric_limits<int>::max();
  if( G.cols() >= mx || G.rows() >= mx ||
      CE.rows() >= mx || CE.cols() >= mx ||
      CI.rows() >= mx || CI.cols() >= mx ||
      ci0.size() >= mx || ce0.size() >= mx || g0.size() >= mx )
    {
    msg << "The dimensions of one of the input matrices or vectors were "
        << "too large." << std::endl
        << "The maximum allowable size for inputs to solve_quadprog is:"
        << mx << std::endl;
    throw std::logic_error(msg.str());
    }

  const int n = static_cast<int>(G.cols());
  const int p = static_cast<int>(CE.cols());
  const int m = static_cast<int>(CI.cols());

  Eigen::LLT<ukfMatrixType, Eigen::Lower> chol(G.cols());

  if( (int)G.rows() != n )
    {
    msg << "The matrix G is not a square matrix (" << G.rows() << " x "
        << G.cols() << ")";
    throw std::logic_error(msg.str() );
    }
  if( (int)CE.rows() != n )
    {
    msg << "The matrix CE is incompatible (incorrect number of rows "
        << CE.rows() << " , expecting " << n << ")";
    throw std::logic_error(msg.str() );
    }
  if( (int)ce0.size() != p )
    {
    msg << "The vector ce0 is incompatible (incorrect dimension "
        << ce0.size() << ", expecting " << p << ")";
    throw std::logic_error(msg.str() );
    }
  if( (int)CI.rows() != n )
    {
    msg << "The matrix CI is incompatible (incorrect number of rows "
        << CI.rows() << " , expecting " << n << ")";
    throw std::logic_error(msg.str() );
    }
  if( (int)ci0.size() != m )
    {
    msg << "The vector ci0 is incompatible (incorrect dimension "
        << ci0.size() << ", expecting " << m << ")";
    throw std::logic_error(msg.str() );
    }
  x.resize(n);
  ukfMatrixType R(n, n), J(n, n);
  ukfVectorType s(m + p), z(n), r(m + p), d(n), np(n), u(m + p), x_old(n), u_old(m + p);

  ukfPrecisionType ss = ukfZero;
  ukfPrecisionType             inf;
  if( std::numeric_limits<ukfPrecisionType>::has_infinity )
    {
    inf = std::numeric_limits<ukfPrecisionType>::infinity();
    }
  else
    {
    inf = 1.0E300;
    }
  Eigen::VectorXi          A(m + p), A_old(m + p), iai(m + p);
  int                      iq, iter = 0;
  Eigen::Matrix<unsigned int, Eigen::Dynamic,1> iaexcl(m + p);

  /* p is the number of equality constraints */
  /* m is the number of inequality constraints */
#ifdef TRACE_SOLVER
  std::cout << std::endl << "Starting solve_quadprog" << std::endl;
  std::cout << "G" << G << std::endl;
  std::cout << "g0" << g0 << std::endl;
  std::cout << "CE" << CE << std::endl;
  std::cout << "ce0" << ce0 << std::endl;
  std::cout << "CI" << CI << std::endl;
  std::cout << "ci0" << ci0 << std::endl;
#endif

  /*
   * Preprocessing phase
   */

  /* compute the trace of the original matrix G */
  ukfPrecisionType c1 = G.trace();

  /* decompose the matrix G in the form L^T L */
  chol.compute(G);

#ifdef TRACE_SOLVER
  std::cout << "G" << G;
#endif
  /* initialize the matrix R */
  d.setZero();
  R.setZero();
  ukfPrecisionType R_norm = ukfOne; /* this variable will hold the norm of the matrix R */

  /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
  J.setIdentity();
  J = chol.matrixU().solve(J);
  ukfPrecisionType c2 = J.trace();
#ifdef TRACE_SOLVER
  std::cout << "J" << J << std::endl;
#endif

  /* c1 * c2 is an estimate for cond(G) */

  /*
   * Find the unconstrained minimizer of the quadratic form ukfHalf * x G x + g0 x
   * this is a feasible point in the dual space
   * x = G^-1 * g0
   */
  x = chol.solve(g0);
  x = -x;
  /* and compute the current solution value */
  ukfPrecisionType f_value = ukfHalf * g0.dot(x);
#ifdef TRACE_SOLVER
  std::cout << "Unconstrained solution: " << f_value << std::endl;
  std::cout << "x" << x << std::endl;
#endif

  /* Add equality constraints to the working set A */
  iq = 0;
  for( int i = 0; i < p; ++i )
    {
    np = CE.col(i);
    compute_d(d, J, np);
    update_z(z, J, d, iq);
    update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
    std::cout << "R" << R << std::endl;
    std::cout << "z" << z << std::endl;
    std::cout << "r" << r << std::endl;
    std::cout << "d" << d << std::endl;
#endif

    /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
     *      becomes feasible */
    ukfPrecisionType t2 = ukfZero;
    if (abs(z.dot(z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
        t2 = (-np.dot(x) - ce0(i)) / z.dot(np);
    
    /* set x = x + t2 * z */
    x += t2 * z;

    /* set u = u+ */
    u(iq) = t2;
    u.head(iq) -= t2 * r.head(iq);

    /* compute the new solution value */
    f_value += ukfHalf * (t2 * t2) * z.dot(np);
    A(i) = -i - 1;

    // NOTE: Removed by CB, its okay not to add an equality constraint!
    // if (!add_constraint(R, J, d, iq, R_norm))
    // {
    // Equality constraints are linearly dependent
    // throw std::runtime_error("Equality constraints are linearly dependent");
    // return f_value;
    // }
    }
  /* set iai = K \ A */
  for( int i = 0; i < m; ++i )
      iai(i) = i;

l1:
  iter++;
#ifdef TRACE_SOLVER
  std::cout << "x" << x << std::endl;
#endif
  /* step 1: choose a violated constraint */
  for( int i = p; i < iq; ++i )
    {
    int ip = A(i);
    iai(ip) = -1;
    }

  /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
  ukfPrecisionType sum = ukfZero;
  ukfPrecisionType psi = ukfZero; /* this value will contain the sum of all infeasibilities */
  int ip = 0; /* ip will be the index of the chosen violated constraint */
  for( int i = 0; i < m; ++i )
    {
      iaexcl[i] = true;
      sum = CI.col(i).dot(x) + ci0(i);
      s(i) = sum;
      psi += std::min(0.0, sum);
    }
#ifdef TRACE_SOLVER
  std::cout << "s" << s << std::endl;
#endif

  if (abs(psi) <= m * std::numeric_limits<double>::epsilon() * c1 * c2 * 100.0)
    {
    /* numerically there are not infeasibilities anymore */
    return f_value;
    }
  /* save old values for u, x and A */
  u_old.head(iq) = u.head(iq);
  A_old.head(iq) = A.head(iq);
  x_old = x;

l2: /* Step 2: check for feasibility and determine a new S-pair */
  for( int i = 0; i < m; ++i )
    {
      if (s(i) < ss && iai(i) != -1 && iaexcl[i])
      {
          ss = s(i);
          ip = i;
      }
    }
  if( ss >= ukfZero )
    {
    return f_value;
    }
  /* set np = n[ip] */
  np = CI.col(ip);
  /* set u = [u 0]^T */
  u(iq) = ukfZero;
  /* add ip to the active set A */
  A(iq) = ip;

#ifdef TRACE_SOLVER
  std::cout << "Trying with constraint " << ip << std::endl;
  std::cout << "np" << np << std::endl;
#endif

l2a: /* Step 2a: determine step direction */
     /* compute z = H np: the step direction in the primal space (through J, see the paper) */
  compute_d(d, J, np);
  update_z(z, J, d, iq);
  /* compute N* np (if q > 0): the negative of the step direction in the dual space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cout << "Step direction z" << std::endl;
  std::cout << "z" << z << std::endl;
  std::cout << "r" << r << std::endl;
  std::cout << "u" << u << std::endl;
  std::cout << "d" << d << std::endl;
  std::cout << "A" << A << std::endl;
#endif

  /* Step 2b: compute step length */
  unsigned int l = 0;
  /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
  ukfPrecisionType t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+(x) / r */
  for( int k = p; k < iq; k++ )
    {
      double tmp;
      if (r(k) > 0.0 && ((tmp = u(k) / r(k)) < t1))
      {
          t1 = tmp;
          l = A(k);
      }
    }
  ukfPrecisionType t2 = inf;
  /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
  if (abs(z.dot(z)) > std::numeric_limits<double>::epsilon()) // i.e. z != 0
      t2 = -s(ip) / z.dot(np);

  /* the step is chosen as the minimum of t1 and t2 */
  ukfPrecisionType t = std::min(t1, t2);
#ifdef TRACE_SOLVER
  std::cout << "Step sizes: " << t << " (t1 = " << t1 << ", t2 = " << t2 << ") ";
#endif

  /* Step 2c: determine new S-pair and take step: */

  /* case (i): no step in primal or dual space */
  if( t >= inf )
    {
    /* QPP is infeasible */
    // FIXME: unbounded to raise
    return inf;
    }
  /* case (ii): step in dual space */
  if( t2 >= inf )
    {
    /* set u = u +  t * [-r 1) and drop constraint l from the active set A */
    u.head(iq) -= t * r.head(iq);
    u(iq) += t;
    iai(l) = l;
    delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
    std::cout << " in dual space: "
              << f_value << std::endl;
    std::cout << "x" << x << std::endl;
    std::cout << "z" << z << std::endl;
    std::cout << "A" << A << std::endl;
#endif
    goto l2a;
    }
  /* case (iii): step in primal and dual space */
  /* set x = x + t * z */
  x += t * z;
  /* update the solution value */
  f_value += t * z.dot(np) * (ukfHalf * t + u(iq));
  /* u = u + t * [-r 1] */
  u.head(iq) -= t * r.head(iq);
  u(iq) += t;
#ifdef TRACE_SOLVER
  std::cout << " in both spaces: "
            << f_value << std::endl;
  std::cout << "x" << x << std::endl;
  std::cout << "u" << u << std::endl;
  std::cout << "r" << r << std::endl;
  std::cout << "A" << A << std::endl;
#endif

  if( fabs(t - t2) < std::numeric_limits<ukfPrecisionType>::epsilon() )
    {
#ifdef TRACE_SOLVER
    std::cout << "Full step has taken " << t << std::endl;
    std::cout << "x" << x << std::endl;
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
    if( !add_constraint(R, J, d, iq, R_norm) )
      {
      iaexcl[ip] = false;
      delete_constraint(R, J, A, u, p, iq, ip);
#ifdef TRACE_SOLVER
    std::cout << "R" << R << std::endl;
    std::cout << "A" << A << std::endl;
#endif
      for( int i = 0; i < m; ++i )
          iai(i) = i;
      for( int i = p; i < iq; ++i )
        {
          A(i) = A_old(i);
          iai(A(i)) = -1;
          u(i) = u_old(i);
        }
      x = x_old;
      goto l2; /* go to step 2 */
      }
    else
        iai(ip) = -1;
#ifdef TRACE_SOLVER
    std::cout << "R" << R << std::endl;
    std::cout << "A" << A << std::endl;
#endif
    goto l1;
    }

  /* a patial step has taken */
#ifdef TRACE_SOLVER
  std::cout << "Partial step has taken " << t << std::endl;
  std::cout << "x" << x << std::endl;
#endif
  /* drop constraint l */
  iai(l) = l;
  delete_constraint(R, J, A, u, p, iq, l);
#ifdef TRACE_SOLVER
  std::cout << "R" << R << std::endl;
  std::cout << "A" << A << std::endl;
#endif

  /* update s[ip] = CI * x + ci0 */
  s(ip) = CI.col(ip).dot(x) + ci0(ip);

#ifdef TRACE_SOLVER
  std::cout << "s" << s << std::endl;
#endif
  goto l2a;
}

}

/*
 *
 * Author: Luca Di Gaspero
 * DIEGM - University of Udine, Italy
 * l.digaspero@uniud.it
 * http://www.diegm.uniud.it/digaspero/
 *
 * LICENSE
 *
 * This file is part of QuadProg++: a C++ library implementing
 * the algorithm of Goldfarb and Idnani for the solution of a (convex)
 * Quadratic Programming problem by means of an active-set dual method.
 * Copyright (C) 2007-2009 Luca Di Gaspero.
 * Copyright (C) 2009 Eric Moyer.
 *
 * QuadProg++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QuadProg++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with QuadProg++. If not, see <http://www.gnu.org/licenses/>.
 *
 */
