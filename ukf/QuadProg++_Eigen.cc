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
  const int    n = static_cast<int>(d.size());

  /* compute d = H^T * np */
  for( int c = 0; c < n; ++c )
    {
    ukfPrecisionType sum = ukfZero;
    for( int r = 0; r < n; ++r )
      {
      sum += J(r,c) * np[r];
      }
    d[c] = sum;
    }
}

inline void update_z(ukfVectorType& z, const ukfMatrixType& J, const ukfVectorType& d, int iq)
{
  const int n = static_cast<int>(z.size());

  /* setting of z = H * d */
  for( int i = 0; i < n; ++i )
    {
    z[i] = ukfZero;
    for( int j = iq; j < n; ++j )
      {
      z[i] += J(i,j) * d[j];
      }
    }
}

inline void update_r(const ukfMatrixType& R, ukfVectorType& r, const ukfVectorType& d, int iq)
{
  /* setting of r = R^-1 d */
  for( int i = iq - 1; i >= 0; i-- )
    {
    ukfPrecisionType sum = ukfZero;
    for( int j = i + 1; j < iq; ++j )
      {
      sum += R(i,j) * r[j];
      }
    r[i] = (d[i] - sum) / R(i,i);
    }
}

inline ukfPrecisionType distance(ukfPrecisionType a, ukfPrecisionType b)
{
  const ukfPrecisionType a1 = fabs(a);
  const ukfPrecisionType b1 = fabs(b);
  if( a1 > b1 )
    {
    ukfPrecisionType t = (b1 / a1);
    return a1 * ::std::sqrt(ukfOne + t * t);
    }
  else if( b1 > a1 )
    {
    ukfPrecisionType t = (a1 / b1);
    return b1 * ::std::sqrt(ukfOne + t * t);
    }
  return a1 * ::std::sqrt(2.0);
}

bool add_constraint(ukfMatrixType& R, ukfMatrixType& J, ukfVectorType& d, int& iq, ukfPrecisionType& R_norm)
{
  const int n = static_cast<int>(d.size());

#ifdef TRACE_SOLVER
  std::cout << "Add constraint " << iq << '/';
#endif
  /* we have to find the Givens rotation which will reduce the element
   *    d[j] to zero.
   *    if it is already zero we don't have to do anything, except of
   *    decreasing j */
  for( int j = n - 1; j >= iq + 1; j-- )
    {
    /* The Givens rotation is done with the matrix (cc cs, cs -cc).
     *    If cc is one, then element (j) of d is zero compared with element
     *    (j - 1). Hence we don't have to do anything.
     *    If cc is zero, then we just have to switch column (j) and column (j - 1)
     *    of J. Since we only switch columns in J, we have to be careful how we
     *    update d depending on the sign of gs.
     *    Otherwise we have to apply the Givens rotation to these columns.
     *    The i - 1 element of d has to be updated to h. */
    ukfPrecisionType cc = d[j - 1];
    ukfPrecisionType ss = d[j];
    ukfPrecisionType h = distance(cc, ss);
    if( fabs(h) < std::numeric_limits<ukfPrecisionType>::epsilon() ) // h == 0
      {
      continue;
      }
    d[j] = ukfZero;
    ss = ss / h;
    cc = cc / h;
    if( cc < ukfZero )
      {
      cc = -cc;
      ss = -ss;
      d[j - 1] = -h;
      }
    else
      {
      d[j - 1] = h;
      }
    ukfPrecisionType xny = ss / (ukfOne + cc);
    for( int k = 0; k < n; k++ )
      {
      ukfPrecisionType t1 = J(k,j - 1);
      ukfPrecisionType t2 = J(k,j);
      J(k,j - 1) = t1 * cc + t2 * ss;
      J(k,j) = xny * (t1 + J(k,j - 1)) - t2;
      }
    }
  /* update the number of constraints added*/
  iq++;
  /* To update R we have to put the iq components of the d vector
   *    into column iq - 1 of R
   */
  for( int i = 0; i < iq; ++i )
    {
    R(i,iq - 1) = d[i];
    }
#ifdef TRACE_SOLVER
  std::cout << iq << std::endl;
  print_matrix("R", R, iq, iq);
  print_matrix("J", J);
  print_vector("d", d, iq);
#endif

  if( fabs(d[iq - 1]) <= std::numeric_limits<ukfPrecisionType>::epsilon() * R_norm )
    {
    // problem degenerate
    return false;
    }
  R_norm = std::max<ukfPrecisionType>(R_norm, fabs(d[iq - 1]) );
  return true;
}

void delete_constraint(ukfMatrixType& R, ukfMatrixType& J, Eigen::VectorXi& A, ukfVectorType& u, int n,
                       int p, int& iq, int l)
{
#ifdef TRACE_SOLVER
  std::cout << "Delete constraint " << l << ' ' << iq;
#endif
  int qq = -1; // just to prevent warnings from smart compilers
  /* Find the index qq for active constraint l to be removed */
  for( int i = p; i < iq; ++i )
    {
    if( A[i] == l )
      {
      qq = i;
      break;
      }
    }
  /* remove the constraint from the active set and the duals */
  for( int i = qq; i < iq - 1; ++i )
    {
    A[i] = A[i + 1];
    u[i] = u[i + 1];
    for( int j = 0; j < n; ++j )
      {
      R(j,i) = R(j,i + 1);
      }
    }

  A[iq - 1] = A[iq];
  u[iq - 1] = u[iq];
  A[iq] = 0;
  u[iq] = ukfZero;
  for( int j = 0; j < iq; ++j )
    {
    R(j,iq - 1) = ukfZero;
    }
  /* constraint has been fully removed */
  iq--;
#ifdef TRACE_SOLVER
  std::cout << '/' << iq << std::endl;
#endif

  if( iq == 0 )
    {
    return;
    }
  for( int j = qq; j < iq; ++j )
    {
    ukfPrecisionType cc = R(j,j);
    ukfPrecisionType ss = R(j + 1,j);
    ukfPrecisionType h = distance(cc, ss);
    if( fabs(h) < std::numeric_limits<ukfPrecisionType>::epsilon() ) // h == 0
      {
      continue;
      }
    cc = cc / h;
    ss = ss / h;
    R(j + 1,j) = ukfZero;
    if( cc < ukfZero )
      {
      R(j,j) = -h;
      cc = -cc;
      ss = -ss;
      }
    else
      {
      R(j,j) = h;
      }

    ukfPrecisionType xny = ss / (ukfOne + cc);
    for( int k =j + 1; k < iq; ++k )
      {
      ukfPrecisionType t1 = R(j,k);
      ukfPrecisionType t2 = R(j + 1,2);
      R(j,k) = t1 * cc + t2 * ss;
      R(j + 1,2) = xny * (t1 + R(j,k)) - t2;
      }
    for( int k = 0; k < n; ++k )
      {
      ukfPrecisionType t1 = J(k,j);
      ukfPrecisionType t2 = J(k,j + 1);
      J(k,j) = t1 * cc + t2 * ss;
      J(k,j + 1) = xny * (J(k,j) + t1) - t2;
      }
    }
}

void cholesky_decomposition(ukfMatrixType& A)
{
  const int n = static_cast<int>(A.rows());

  for( int i = 0; i < n; ++i )
    {
    for( int j = i; j < n; ++j )
      {
      ukfPrecisionType sum = A(i,j);
      for( int k =i - 1; k >= 0; k-- )
        {
        sum -= A(i,k) * A(j,k);
        }
      if( i == j )
        {
        if( sum <= ukfZero )
          {
          std::ostringstream os;
          // raise error
          std::cout << "A" << A;
          os << "Error in cholesky decomposition, sum: " << sum;
          throw std::logic_error(os.str());
          }
        A(i,i) = ::std::sqrt(sum);
        }
      else
        {
        A(j,i) = sum / A(i,i);
        }
      }
    for( int k =i + 1; k < n; ++k )
      {
      A(i,k) = A(k,i);
      }
    }
}

inline void forward_elimination(const ukfMatrixType& L, ukfVectorType& y, const ukfVectorType& b)
{
  const int n = static_cast<int>(L.rows());

  y[0] = b[0] / L(0,0);
  for( int i = 1; i < n; ++i )
    {
    y[i] = b[i];
    for( int j = 0; j < i; ++j )
      {
      y[i] -= L(i,j) * y[j];
      }
    y[i] = y[i] / L(i,i);
    }
}

inline void backward_elimination(const ukfMatrixType& U, ukfVectorType& x, const ukfVectorType& y)
{
  const int n = static_cast<int>(U.rows());

  x[n - 1] = y[n - 1] / U(n - 1,n - 1);
  for( int i = n - 2; i >= 0; i-- )
    {
    x[i] = y[i];
    for( int j = i + 1; j < n; ++j )
      {
      x[i] -= U(i,j) * x[j];
      }
    x[i] = x[i] / U(i,i);
    }
}

void cholesky_solve(const ukfMatrixType& L, ukfVectorType& x, const ukfVectorType& b)
{
  const int n = static_cast<int>(L.rows());

  ukfVectorType y(n);

  /* Solve L * y = b */
  forward_elimination(L, y, b);
  /* Solve L^T * x = y */
  backward_elimination(L, x, y);
}


inline ukfPrecisionType dot_product(const ukfVectorType& x, const ukfVectorType& y)
{
  int    n = static_cast<int>(x.size());
  ukfPrecisionType sum;

  sum = ukfZero;
  for( int i = 0; i < n; ++i )
    {
    sum += x[i] * y[i];
    }
  return sum;
}

// Utility functions for printing vectors and matrices
void print_matrix(const char* name, const ukfMatrixType& A, int n, int m)
{
  std::ostringstream s;
  std::string        t;

  if( n == -1 )
    {
    n = static_cast<int>(A.rows());
    }
  if( m == -1 )
    {
    m = static_cast<int>(A.cols());
    }

  s << name << ": " << std::endl;
  for( int i = 0; i < n; ++i )
    {
    s << " ";
    for( int j = 0; j < m; ++j )
      {
      s << A(i,j) << ", ";
      }
    s << std::endl;
    }
  t = s.str();
  t = t.substr(0, t.size() - 3); // To remove the trailing space, comma and newline

  std::cout << t << std::endl;
}

template <typename T>
void print_vector(const char* name, const typename Eigen::Matrix<T,Eigen::Dynamic, 1> & v, int n)
{
  std::ostringstream s;
  std::string        t;

  if( n == -1 )
    {
    n = v.size();
    }

  s << name << ": " << std::endl << " ";
  for( int i = 0; i < n; ++i )
    {
    s << v[i] << ", ";
    }
  t = s.str();
  t = t.substr(0, t.size() - 2); // To remove the trailing space and comma

  std::cout << t << std::endl;
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
  print_matrix("G", G);
  print_vector("g0", g0);
  print_matrix("CE", CE);
  print_vector("ce0", ce0);
  print_matrix("CI", CI);
  print_vector("ci0", ci0);
#endif

  /*
   * Preprocessing phase
   */

  /* compute the trace of the original matrix G */
  ukfPrecisionType c1 = ukfZero;
  for(int i = 0; i < n; ++i )
    {
    c1 += G(i,i);
    }
  /* decompose the matrix G in the form L^T L */
  cholesky_decomposition(G);
#ifdef TRACE_SOLVER
  print_matrix("G", G);
#endif
  /* initialize the matrix R */
  for( int i = 0; i < n; ++i )
    {
    d[i] = ukfZero;
    for( int j = 0; j < n; ++j )
      {
      R(i,j) = ukfZero;
      }
    }
  ukfPrecisionType R_norm = ukfOne; /* this variable will hold the norm of the matrix R */

  /* compute the inverse of the factorized matrix G^-1, this is the initial value for H */
  ukfPrecisionType c2 = ukfZero;
  for( int i = 0; i < n; ++i )
    {
    d[i] = ukfOne;
    forward_elimination(G, z, d);
    for( int j = 0; j < n; ++j )
      {
      J(i,j) = z[j];
      }
    c2 += z[i];
    d[i] = ukfZero;
    }
#ifdef TRACE_SOLVER
  print_matrix("J", J);
#endif

  /* c1 * c2 is an estimate for cond(G) */

  /*
   * Find the unconstrained minimizer of the quadratic form ukfHalf * x G x + g0 x
   * this is a feasible point in the dual space
   * x = G^-1 * g0
   */
  cholesky_solve(G, x, g0);
  for( int i = 0; i < n; ++i )
    {
    x[i] = -x[i];
    }
  /* and compute the current solution value */
  ukfPrecisionType f_value = ukfHalf * dot_product(g0, x);
#ifdef TRACE_SOLVER
  std::cout << "Unconstrained solution: " << f_value << std::endl;
  print_vector("x", x);
#endif

  /* Add equality constraints to the working set A */
  iq = 0;
  for( int i = 0; i < p; ++i )
    {
    for( int j = 0; j < n; ++j )
      {
      np[j] = CE(j,i);
      }
    compute_d(d, J, np);
    update_z(z, J, d, iq);
    update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
    print_matrix("R", R, n, iq);
    print_vector("z", z);
    print_vector("r", r, iq);
    print_vector("d", d);
#endif

    /* compute full step length t2: i.e., the minimum step in primal space s.t. the contraint
     *      becomes feasible */
    ukfPrecisionType t2 = ukfZero;
    if( fabs(dot_product(z, z) ) > std::numeric_limits<ukfPrecisionType>::epsilon() ) // i.e. z != 0
      {
      t2 = (-dot_product(np, x) - ce0[i]) / dot_product(z, np);
      }
    /* set x = x + t2 * z */
    for( int k = 0; k < n; ++k )
      {
      x[k] += t2 * z[k];
      }

    /* set u = u+ */
    u[iq] = t2;
    for( int k = 0; k < iq; ++k )
      {
      u[k] -= t2 * r[k];
      }

    /* compute the new solution value */
    f_value += ukfHalf * (t2 * t2) * dot_product(z, np);
    A[i] = -i - 1;

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
    {
    iai[i] = i;
    }

l1:
  iter++;
#ifdef TRACE_SOLVER
  print_vector("x", x);
#endif
  /* step 1: choose a violated constraint */
  for( int i = p; i < iq; ++i )
    {
    int ip = A[i];
    iai[ip] = -1;
    }

  /* compute s[x] = ci^T * x + ci0 for all elements of K \ A */
  ukfPrecisionType ss = ukfZero;
  ukfPrecisionType psi = ukfZero; /* this value will contain the sum of all infeasibilities */
  int ip = 0; /* ip will be the index of the chosen violated constraint */
  for( int i = 0; i < m; ++i )
    {
    iaexcl[i] = true;
    ukfPrecisionType sum = ukfZero;
    for( int j = 0; j < n; ++j )
      {
      sum += CI(j,i) * x[j];
      }
    sum += ci0[i];
    s[i] = sum;
    psi += std::min(ukfZero, sum);
    }
#ifdef TRACE_SOLVER
  print_vector("s", s, m);
#endif

  if( fabs(psi) <= m * std::numeric_limits<ukfPrecisionType>::epsilon() * c1 * c2 * 100.0 )
    {
    /* numerically there are not infeasibilities anymore */
    return f_value;
    }
  /* save old values for u and A */
  for( int i = 0; i < iq; ++i )
    {
    u_old[i] = u[i];
    A_old[i] = A[i];
    }
  /* and for x */
  for( int i = 0; i < n; ++i )
    {
    x_old[i] = x[i];
    }

l2: /* Step 2: check for feasibility and determine a new S-pair */
  for( int i = 0; i < m; ++i )
    {
    if( s[i] < ss && iai[i] != -1 && iaexcl[i] )
      {
      ss = s[i];
      ip = i;
      }
    }
  if( ss >= ukfZero )
    {
    return f_value;
    }
  /* set np = n[ip] */
  for( int i = 0; i < n; ++i )
    {
    np[i] = CI(i,ip);
    }
  /* set u = [u 0]^T */
  u[iq] = ukfZero;
  /* add ip to the active set A */
  A[iq] = ip;

#ifdef TRACE_SOLVER
  std::cout << "Trying with constraint " << ip << std::endl;
  print_vector("np", np);
#endif

l2a: /* Step 2a: determine step direction */
     /* compute z = H np: the step direction in the primal space (through J, see the paper) */
  compute_d(d, J, np);
  update_z(z, J, d, iq);
  /* compute N* np (if q > 0): the negative of the step direction in the dual space */
  update_r(R, r, d, iq);
#ifdef TRACE_SOLVER
  std::cout << "Step direction z" << std::endl;
  print_vector("z", z);
  print_vector("r", r, iq + 1);
  print_vector("u", u, iq + 1);
  print_vector("d", d);
  print_vector("A", A, iq + 1);
#endif

  /* Step 2b: compute step length */
  unsigned int l = 0;
  /* Compute t1: partial step length (maximum step in dual space without violating dual feasibility */
  ukfPrecisionType t1 = inf; /* +inf */
  /* find the index l s.t. it reaches the minimum of u+[x] / r */
  for( int k = p; k < iq; k++ )
    {
    if( r[k] > ukfZero )
      {
      if( u[k] / r[k] < t1 )
        {
        t1 = u[k] / r[k];
        l = A[k];
        }
      }
    }
  ukfPrecisionType t2 = inf;
  /* Compute t2: full step length (minimum step in primal space such that the constraint ip becomes feasible */
  if( fabs(dot_product(z, z) )  > std::numeric_limits<ukfPrecisionType>::epsilon() ) // i.e. z != 0
    {
    t2 = -s[ip] / dot_product(z, np);
    }

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
    /* set u = u +  t * [-r 1] and drop constraint l from the active set A */
    for( int k = 0; k < iq; k++ )
      {
      u[k] -= t * r[k];
      }
    u[iq] += t;
    iai[l] = l;
    delete_constraint(R, J, A, u, n, p, iq, l);
#ifdef TRACE_SOLVER
    std::cout << " in dual space: "
              << f_value << std::endl;
    print_vector("x", x);
    print_vector("z", z);
    print_vector("A", A, iq + 1);
#endif
    goto l2a;
    }
  /* case (iii): step in primal and dual space */
  /* set x = x + t * z */
  for( int k = 0; k < n; k++ )
    {
    x[k] += t * z[k];
    }
  /* update the solution value */
  f_value += t * dot_product(z, np) * (ukfHalf * t + u[iq]);
  /* u = u + t * [-r 1] */
  for( int k = 0; k < iq; k++ )
    {
    u[k] -= t * r[k];
    }
  u[iq] += t;
#ifdef TRACE_SOLVER
  std::cout << " in both spaces: "
            << f_value << std::endl;
  print_vector("x", x);
  print_vector("u", u, iq + 1);
  print_vector("r", r, iq + 1);
  print_vector("A", A, iq + 1);
#endif

  if( fabs(t - t2) < std::numeric_limits<ukfPrecisionType>::epsilon() )
    {
#ifdef TRACE_SOLVER
    std::cout << "Full step has taken " << t << std::endl;
    print_vector("x", x);
#endif
    /* full step has taken */
    /* add constraint ip to the active set*/
    if( !add_constraint(R, J, d, iq, R_norm) )
      {
      iaexcl[ip] = false;
      delete_constraint(R, J, A, u, n, p, iq, ip);
#ifdef TRACE_SOLVER
      print_matrix("R", R);
      print_vector("A", A, iq);
      print_vector("iai", iai);
#endif
      for( int i = 0; i < m; ++i )
        {
        iai[i] = i;
        }
      for( int i = p; i < iq; ++i )
        {
        A[i] = A_old[i];
        u[i] = u_old[i];
        iai[A[i]] = -1;
        }
      for( int i = 0; i < n; ++i )
        {
        x[i] = x_old[i];
        }
      goto l2; /* go to step 2 */
      }
    else
      {
      iai[ip] = -1;
      }
#ifdef TRACE_SOLVER
    print_matrix("R", R);
    print_vector("A", A, iq);
    print_vector("iai", iai);
#endif
    goto l1;
    }

  /* a patial step has taken */
#ifdef TRACE_SOLVER
  std::cout << "Partial step has taken " << t << std::endl;
  print_vector("x", x);
#endif
  /* drop constraint l */
  iai[l] = l;
  delete_constraint(R, J, A, u, n, p, iq, l);
#ifdef TRACE_SOLVER
  print_matrix("R", R);
  print_vector("A", A, iq);
#endif

  /* update s[ip] = CI * x + ci0 */
  ukfPrecisionType sum = ukfZero;
  for( int k = 0; k < n; k++ )
    {
    sum += CI(k,ip) * x[k];
    }
  s[ip] = sum + ci0[ip];

#ifdef TRACE_SOLVER
  print_vector("s", s, m);
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
