#include "BPDN_homotopy.h"
#include <vector>
#include <limits>

template < typename T>
T sign(T val)
{
  if(val > 0)
    {
    return 1.0;
    }
  else if(val < 0)
    {
    return -1.0;
    }
  return 0.0;
}
typedef std::vector<unsigned long> IndexVec;

// function [i_delta, out_x, delta, chk_x] = update_primal(gamma_x,
// gamma_lambda, z_x, x_k, del_x_vec, pk, dk, ep// silon, out_lambda);
void
update_primal(IndexVec & gamma_x, // current support of x
              IndexVec & gamma_lambda, // current support of
                                                         // lambda
              MatrixType &                 z_x, // sign sequence of x
              MatrixType &                 x_k, // sign sequence of lambda
              MatrixType &                 del_x_vec, // primal update direction
              MatrixType &                 pk,       //
              MatrixType &                 dk,       //
              double                       epsilon, // current value of epsilon
              IndexVec & out_lambda, // element removed from
              // support of lambda in
              // previous step if any
              // OUTPUTS
              unsigned long &              i_delta, // index corresponding to
              // newly active primal
              // constraint (new_lambda)
              IndexVec & out_x, // element in x shrunk to zero;
              double  &                    delta, // primal step size
              unsigned &                   chk_x) // 1 an element is removed
                                                  // from support of x
                                                  // 0 a new element enters
                                                  // the support of lambd
{
//
// N = length(x_k);
  const unsigned long N(x_k.rows());
//
// temp_gamma = zeros(N,1);
  IndexVec temp_gamma(N,0);
// temp_gamma(gamma_lambda) = gamma_lambda;
  for(unsigned int i = 0; i < gamma_lambda.size(); ++i)
    {
    temp_gamma[gamma_lambda[i]] = gamma_lambda[i];
    }
// gamma_lc = find([1:N]' ~= temp_gamma);
  IndexVec gamma_lc;
  for(unsigned int i = 0; i < temp_gamma.size(); ++i)
    {
    if(temp_gamma[i] == 0)
      {
      gamma_lc.push_back(i);
      }
    }
//
// delta1_constr = (epsilon-pk(gamma_lc))./(1+dk(gamma_lc));
  MatrixType delta1_constr(gamma_lc.size(),1);
  for(unsigned i = 0; i < gamma_lc.size(); ++i)
    {
    delta1_constr(i,0) = (epsilon - pk(gamma_lc[i]))/(1.0 + dk(gamma_lc[i]));
    }
// delta1_pos_ind = find(delta1_constr>0);
  IndexVec delta1_pos_ind;
  for(unsigned int i = 0; i < delta1_constr.rows(); ++i)
    {
    if(delta1_constr(i,0) > 0.0)
      {
      delta1_pos_ind.push_back(i);
      }
    }
// delta1_pos = delta1_constr(delta1_pos_ind);
// [delta1 i_delta1] = min(delta1_pos);
  MatrixType delta1_pos(delta1_pos_ind.size(),1);
  double delta1 = std::numeric_limits<double>::max();
  long i_delta1 = -1;
  for(unsigned int i = 0; i < delta1_pos_ind.size(); ++i)
    {
    double cur = delta1_constr(delta1_pos_ind[i],0);
    delta1_pos(i,0) = cur;
    if(cur < delta1)
      {
      delta1 = cur;
      i_delta1 = i;
      }
    }
// if isempty(delta1)
//     delta1 = inf;
// end
  if(delta1 == std::numeric_limits<double>::max())
    {
    delta1 = std::numeric_limits<double>::infinity();
    }
// delta2_constr = (epsilon+pk(gamma_lc))./(1-dk(gamma_lc));
  MatrixType delta2_constr(gamma_lc.size(),1);
  for(unsigned int i = 0; i < gamma_lc.size(); ++i)
    {
    delta2_constr(i,0) = (epsilon + pk(gamma_lc[i],0))/(1.0 - dk(gamma_lc[i]));
    }
// delta2_pos_ind = find(delta2_constr>0);
  IndexVec delta2_pos_ind;
  for(unsigned int i = 0; i < delta2_constr.rows(); ++i)
    {
    if(delta2_constr(i,0) > 0.0)
      {
      delta2_pos_ind.push_back(i);
      }
    }
// delta2_pos = delta2_constr(delta2_pos_ind);
// [delta2 i_delta2] = min(delta2_pos);
  MatrixType delta2_pos(delta2_pos_ind.size(),1);
  double delta2 = std::numeric_limits<double>::max();
  long i_delta2 = -1;
  for(unsigned int i = 0; i < delta2_pos_ind.size(); ++i)
    {
    double cur = delta2_constr(delta2_pos_ind[i],0);
    delta2_pos(i,0) = cur;
    if(cur < delta2)
      {
      delta2 = cur;
      i_delta2 = i;
      }
    }
// if isempty(delta2)
//     delta2 = inf;
// end
  if(delta2 == std::numeric_limits<double>::max())
    {
    delta2 = std::numeric_limits<double>::infinity();
    }

// if delta1>delta2
  if(delta1 > delta2)
    {
//     delta = delta2;
    delta = delta2;
//     i_delta = gamma_lc(delta2_pos_ind(i_delta2));
    i_delta = gamma_lc[delta2_pos_ind[i_delta2]];
    }
// else
  else if(i_delta1 != -1 && delta1_pos_ind.size() > 0) // if delta1_pos_ind is empty
    {
//     delta = delta1;
    delta = delta1;
//     i_delta = gamma_lc(delta1_pos_ind(i_delta1));
    i_delta = gamma_lc[delta1_pos_ind[i_delta1]];
// end
    }

// delta3_constr = (-x_k(gamma_x)./del_x_vec(gamma_x));
  MatrixType delta3_constr(gamma_x.size(),1);
  for(unsigned int i = 0; i < gamma_x.size(); ++i)
    {
    delta3_constr(i,0) = (-x_k(gamma_x[i],0) / del_x_vec(gamma_x[i]));
    }
// delta3_pos_index = find(delta3_constr>0);
  IndexVec delta3_pos_index;
  for(unsigned i = 0; i < delta3_constr.rows(); ++i)
    {
    if(delta3_constr(i) > 0)
      {
      delta3_pos_index.push_back(i);
      }
    }
// [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
  double delta3 = std::numeric_limits<double>::max();
  long i_delta3(-1);
  for(unsigned int i = 0; i < delta3_pos_index.size(); ++i)
    {
    const double cur = delta3_constr(delta3_pos_index[i],0);
    if(cur < delta3)
      {
      delta3 = cur;
      i_delta3 = i;
      }
    }
// out_x_index = gamma_x(delta3_pos_index(i_delta3));
//
// chk_x = 0;
  chk_x = 0;
// out_x = [];
  out_x.clear();
// if delta3 > 0 & delta3 <= delta
  if(delta3 > 0.0 && delta3 <= delta)
//     chk_x = 1;
    {
//     delta = delta3;
    delta = delta3;
    unsigned long out_x_index = gamma_x[delta3_pos_index[i_delta3]];
//     out_x = out_x_index;
    out_x.push_back(out_x_index);
// end
    }
//
// %%% THESE ARE PROBABLY UNNECESSARY
// %%% NEED TO REMOVE THEM.
//
// % The following checks are just to deal with degenerate cases when more
// % than one elements want to enter or leave the support at any step
// % (e.g., Bernoulli matrix with small number of measurements)
//
// % This one is ONLY for those indices which are zero. And we don't know where
// % will its dx point in next steps, so after we calculate dx and its in opposite
// % direction to z_x, we will have to remove that index from the support.
// xk_1 = x_k+delta*del_x_vec;
  MatrixType xk_1 = x_k + delta * del_x_vec;
// xk_1(out_x) = 0;
  for(unsigned int i = 0; i < out_x.size(); ++i)
    {
    xk_1(out_x[i]) = 0.0;
    }
// wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
  IndexVec wrong_sign;
  for(unsigned int i = 0; i < gamma_x.size(); ++i)
    {
    if(sign(xk_1(gamma_x[i],0)) * z_x(gamma_x[i]) == -1)
      {
      wrong_sign.push_back(i);
      }
    }
// if ~isempty(gamma_x(wrong_sign))
  if(wrong_sign.size() > 0)
    {
//     chk_x = 1;
    chk_x = 1;
//     delta = 0;
    delta = 0;
//     % can also choose specific element which became non-zero first but all
//     % that matters here is AtA(gx,gl) doesn't become singular.
//     % [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
//     out_x = gamma_x(wrong_sign(1));
    out_x.clear();
    out_x.push_back(gamma_x[wrong_sign[0]]);
// end
    }
//
// % If more than one primal constraints became active in previous iteration i.e.,
// % more than one elements wanted to enter the support and we added only one.
// % So here we need to check if those remaining elements are still active.
// i_delta_temp = gamma_lc(find(abs(pk(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps));
  IndexVec i_delta_temp;
  for(unsigned int i = 0; i < gamma_lc.size(); ++i)
    {
    if(std::abs(pk(gamma_lc[i],0) + delta * dk(gamma_lc[i],0)) - (epsilon - delta) >= 10 * Eps)
      {
      i_delta_temp.push_back(gamma_lc[i]);
      }
    }
// if ~isempty(i_delta_temp)
  IndexVec i_delta_more;
  if(i_delta_temp.size() > 0)
    {
//     if ~isempty(out_lambda)
    if(out_lambda.size() > 0)
      {
//         i_delta_more = i_delta_temp;%(find(i_delta_temp~=out_lambda));
      i_delta_more = i_delta_temp;
      }
//     else
    else
      {
//         i_delta_more = i_delta_temp;
      i_delta_more = i_delta_temp;
//     end
      }
//     if length(i_delta_more)>=1 & ~sum((i_delta_temp==i_delta))
    unsigned equalcount(0);
    for(unsigned i = 0; i < i_delta_temp.size(); ++i)
      {
      if(i_delta_temp[i] == i_delta)
        {
        ++equalcount;
        }
      }
    if(i_delta_more.size() > 0 && equalcount == 0)
      {
//         % ideal way would be to check that incoming element doesn't make AtA
//         % singular!
//         [v_temp i_temp] = max(-pk(i_delta_more)./dk(i_delta_more));
      double v_temp = std::numeric_limits<double>::min();
      long i_temp = -1;
      for(unsigned i = 0; i < i_delta_more.size(); ++i)
        {
        double curval = -pk(i_delta_more[i],0) / dk(i_delta_more[i],0);
        if(curval > v_temp)
          {
          v_temp = curval;
          i_temp = i;
          }
        }
//         i_delta = i_delta_more(i_temp);
      i_delta = i_delta_more[i_temp];
//         delta = 0;
      delta = 0.0;
//         chk_x = 0;
      chk_x = 0.0;
//         out_x = [];
      out_x.clear();
//     end
      }
// end
    }
}

// % update_primal.m
// %
// % This function computes the minimum step size in the primal update direction and
// % finds change in the primal or dual support with that step.
// %
// % Inputs:
// % gamma_x - current support of x
// % gamma_lambda - current support of lambda
// % z_x - sign sequence of x
// % z_lambda - sign sequence of lambda
// % del_x_vec - primal update direction
// % pk
// % dk
// % epsilon - current value of epsilon
// % out_lambda - element removed from support of lambda in previous step (if any)
// %
// % Outputs:
// % i_delta - index corresponding to newly active primal constraint (new_lambda)
// % out_x - element in x shrunk to zero
// % delta - primal step size
// % chk_x - 1  an element is removed from support of x
// %         0  a new element enters the support of lambda
// %
// % Written by: Salman Asif, Georgia Tech
// % Email: sasif@ece.gatech.edu
//
// function [i_delta, out_x, delta, chk_x] = update_primal(gamma_x,
// gamma_lambda, z_x, x_k, del_x_vec, pk, dk, ep// silon, out_lambda);
//
// N = length(x_k);
//
// % gamma_lc = setdiff([1:N]', [gamma_lambda; out_lambda]); WRONG
// % check out_lambda as well, that is if outgoing lambda switches sign in just one step
// temp_gamma = zeros(N,1);
// temp_gamma(gamma_lambda) = gamma_lambda;
// gamma_lc = find([1:N]' ~= temp_gamma);
//
// delta1_constr = (epsilon-pk(gamma_lc))./(1+dk(gamma_lc));
// delta1_pos_ind = find(delta1_constr>0);
// delta1_pos = delta1_constr(delta1_pos_ind);
// [delta1 i_delta1] = min(delta1_pos);
// if isempty(delta1)
//     delta1 = inf;
// end
// delta2_constr = (epsilon+pk(gamma_lc))./(1-dk(gamma_lc));
// delta2_pos_ind = find(delta2_constr>0);
// delta2_pos = delta2_constr(delta2_pos_ind);
// [delta2 i_delta2] = min(delta2_pos);
// if isempty(delta2)
//     delta2 = inf;
// end

// if delta1>delta2
//     delta = delta2;
//     i_delta = gamma_lc(delta2_pos_ind(i_delta2));
// else
//     delta = delta1;
//     i_delta = gamma_lc(delta1_pos_ind(i_delta1));
// end

// delta3_constr = (-x_k(gamma_x)./del_x_vec(gamma_x));
// delta3_pos_index = find(delta3_constr>0);
// [delta3 i_delta3] = min(delta3_constr(delta3_pos_index));
// out_x_index = gamma_x(delta3_pos_index(i_delta3));
//
// chk_x = 0;
// out_x = [];
// if delta3 > 0 & delta3 <= delta
//     chk_x = 1;
//     delta = delta3;
//     out_x = out_x_index;
// end
//
// %%% THESE ARE PROBABLY UNNECESSARY
// %%% NEED TO REMOVE THEM.
//
// % The following checks are just to deal with degenerate cases when more
// % than one elements want to enter or leave the support at any step
// % (e.g., Bernoulli matrix with small number of measurements)
//
// % This one is ONLY for those indices which are zero. And we don't know where
// % will its dx point in next steps, so after we calculate dx and its in opposite
// % direction to z_x, we will have to remove that index from the support.
// xk_1 = x_k+delta*del_x_vec;
// xk_1(out_x) = 0;
// wrong_sign = find(sign(xk_1(gamma_x)).*z_x(gamma_x)==-1);
// if ~isempty(gamma_x(wrong_sign))
//     chk_x = 1;
//     delta = 0;
//     % can also choose specific element which became non-zero first but all
//     % that matters here is AtA(gx,gl) doesn't become singular.
//     % [val_wrong_x ind_wrong_x] =  sort(abs(del_x_vec(gamma_x(wrong_sign))),'descend');
//     out_x = gamma_x(wrong_sign(1));
// end
//
// % If more than one primal constraints became active in previous iteration i.e.,
// % more than one elements wanted to enter the support and we added only one.
// % So here we need to check if those remaining elements are still active.
// i_delta_temp = gamma_lc(find(abs(pk(gamma_lc)+delta*dk(gamma_lc))-(epsilon-delta) >= 10*eps));
// if ~isempty(i_delta_temp)
//     if ~isempty(out_lambda)
//         i_delta_more = i_delta_temp;%(find(i_delta_temp~=out_lambda));
//     else
//         i_delta_more = i_delta_temp;
//     end
//     if length(i_delta_more)>=1 & ~sum((i_delta_temp==i_delta))
//         % ideal way would be to check that incoming element doesn't make AtA
//         % singular!
//         [v_temp i_temp] = max(-pk(i_delta_more)./dk(i_delta_more));
//         i_delta = i_delta_more(i_temp);
//         delta = 0;
//         chk_x = 0;
//         out_x = [];
//     end
// end

// function [x_out, gamma_x, total_iter, total_time] = BPDN_homotopy_function(A, y, tau, maxiter)
MatrixType BPDN_HOMOTOPY_function(MatrixType &A,MatrixType &y, double tau, unsigned int maxiter)
{
//
// t0 = cputime;
//
// N = size(A,2);
  unsigned long N = A.cols();
// K = size(A,1);
  unsigned long K = A.rows();
//
// % Initialization of primal and dual sign and support
// z_x = zeros(N,1);
  MatrixType z_x;
  z_x = MatrixType::Zero(N,1);
// gamma_x = [];       % Primal support
  IndexVec gamma_x;
//
// % Initial step
// Primal_constrk = -A'*y;
  MatrixType Primal_constrk = (-1 * A.transpose()) * y;
// [c i] = max(abs(Primal_constrk));
  unsigned long irow, icol;
  double c = Primal_constrk.cwiseAbs().maxCoeff(&irow,&icol);
//
// gamma_xk = i;
  IndexVec gamma_xk;
  gamma_xk.push_back(irow);
//
// epsilon = c;
  double epsilon(c);
// xk_1 = zeros(N,1);
  MatrixType xk_1;
  xk_1 = MatrixType::Zero(N,1);
//
// z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
// Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;
  for(unsigned i = 0; i < gamma_xk.size(); ++i)
    {
    z_x(gamma_xk[i],0) = -sign(Primal_constrk(gamma_xk[i],0));
    Primal_constrk(gamma_xk[i],0) = sign(Primal_constrk(gamma_xk[i],0)) * epsilon;
    }

//
// z_xk = z_x;
  MatrixType z_xk = z_x;
//
// % loop parameters
// done = 0;
// iter = 0;
  unsigned long iter = 0;
// data_precision = eps;   % floating point precision
  double data_precision(Eps);
//
// old_delta = 0;
  double old_delta(0.0);
// out_x = [];
  IndexVec out_x;
// count_delta_stop = 0;
  unsigned count_delta_stop = 0;
//
// constraint_plots = 1;
  unsigned constraint_plots = 1;
//
// AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
  MatrixType AtgxAgx;
  {
  MatrixType Atmp(A.rows(),gamma_xk.size());
  for(unsigned int i = 0; i < gamma_xk.size(); ++i)
    {
    Atmp.col(i) = A.col(gamma_xk[i]);
    }
  AtgxAgx = Atmp.transpose() * Atmp;
  }

// iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
  MatrixType iAtgxAgx = AtgxAgx.inverse();
//
// while iter < maxiter
  while(iter < maxiter)
    {
//     iter = iter+1;
    ++iter;
//     % warning('off','MATLAB:divideByZero')
//
//     gamma_x = gamma_xk;
    gamma_x = gamma_xk;
//     z_x = z_xk;
    z_x = z_xk;
//     x_k = xk_1;
    MatrixType x_k = xk_1;
//
//     %%%%%%%%%%%%%%%%%%%%%
//     %%%% update on x %%%%
//     %%%%%%%%%%%%%%%%%%%%%
//
//     % Update direction
//     %del_x = inv(A(:,gamma_x)'*A(:,gamma_x))*z_x(gamma_x);
//     del_x = iAtgxAgx*z_x(gamma_x);
    MatrixType del_x;
    {
    MatrixType ztmp(gamma_x.size(),1);
    for(unsigned i = 0; i < gamma_x.size(); ++i)
      {
      ztmp(i,0) = z_x(gamma_x[i],0);
      }
    del_x = iAtgxAgx * ztmp;
    }
//     del_x_vec = zeros(N,1);
    MatrixType del_x_vec;
    del_x_vec = MatrixType::Zero(N,1);
//     del_x_vec(gamma_x) = del_x;
    for(unsigned i = 0; i < gamma_x.size(); ++i)
      {
      del_x_vec(gamma_x[i],0) = del_x(i,0);
      }
//
//     pk = Primal_constrk;
    MatrixType pk = Primal_constrk;
//     %dk = A'*(A*del_x_vec);
//     Agdelx = A(:,gamma_x)*del_x;
    MatrixType Agdelx;
    {
    MatrixType Atmp(A.rows(),gamma_x.size());
    for(unsigned int i = 0; i < gamma_x.size(); ++i)
      {
      Atmp.col(i) = A.col(gamma_x[i]);
      }
    Agdelx = Atmp * del_x;
    }
//     dk = A'*Agdelx;
    MatrixType dk = A.transpose() * Agdelx;
//
//     %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
//     pk_temp = Primal_constrk;
    MatrixType pk_temp = Primal_constrk;
//     gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
    IndexVec gammaL_temp;
    for(unsigned i = 0; i < Primal_constrk.rows(); ++i)
      {
      if(std::fabs(std::fabs(Primal_constrk(i,0))-epsilon) < std::min(epsilon,2*Eps))
        {
        gammaL_temp.push_back(i);
        }
      }
//     pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
    for(unsigned int i = 0; i < gammaL_temp.size(); ++i)
      {
      pk_temp(gammaL_temp[i],0) = sign(Primal_constrk(gammaL_temp[i],0)) * epsilon;
      }
//
//     xk_temp = x_k;
    MatrixType xk_temp = x_k;
//     gammaX_temp = find(abs(x_k)<2*eps);
//     xk_temp(gammaX_temp) = 0;
    for(unsigned int i = 0; i < x_k.size(); ++i)
      {
      if(std::fabs(x_k(i,0)) < 2*Eps)
        {
        xk_temp(i,0) = 0.0;
        }
      }
//     %%%---
//
//     % Compute the step size
//     [i_delta, out_x, delta, chk_x] = update_primal(gamma_x, gamma_x, z_x,  xk_temp, del_x_vec, pk_temp, dk, epsilon, out_x);
    unsigned long i_delta;
    double        delta;
    unsigned      chk_x;
    update_primal(gamma_x,
                   gamma_x,
                   z_x,
                   xk_temp,
                   del_x_vec,
                   pk_temp,
                   dk,
                   epsilon,
                   out_x,
                   i_delta,
                   out_x,
                   delta,
                   chk_x);
//
//     if old_delta < 4*eps && delta < 4*eps
    if(old_delta < (4.0 * Eps) && delta < (4 * Eps))
      {
//         count_delta_stop = count_delta_stop + 1;
      count_delta_stop++;
      }
//     else
    else
      {
//         count_delta_stop = 0;
      count_delta_stop = 0;
//     end
      }
//     if count_delta_stop >= 500
    if(count_delta_stop >= 500)
      {
//         %disp('stuck in some corner');
//         break;
      break;
//     end
      }
//     old_delta = delta;
    old_delta = delta;
//
//     xk_1 = x_k+delta*del_x_vec;
    xk_1 = x_k + ( delta * del_x_vec);
//     Primal_constrk = pk+delta*dk;
    Primal_constrk = pk + (delta * dk);
//     epsilon_old = epsilon;
    double epsilon_old = epsilon;
//     epsilon = epsilon-delta;
    epsilon -= delta;
//
//     if epsilon <= tau;
    if(epsilon <= tau)
      {
//         xk_1 = x_k + (epsilon_old-tau)*del_x_vec;
      xk_1 = x_k + ((epsilon_old - tau) * del_x_vec);
//         total_time= cputime-t0;
//         break;
      break;
//     end
      }
//
    if(chk_x)
//     if chk_x == 1
      {
//         % If an element is removed from gamma_x
//         gx_old = gamma_x;
      IndexVec gx_old = gamma_x;
//         len_gamma = length(gamma_x);
      unsigned long len_gamma = gamma_x.size();
//
//         outx_index = find(gamma_x==out_x);
      unsigned long outx_index(0);
      for(unsigned int i = 0; i < gamma_x.size(); ++i)
        {
        if(out_x[i] == gamma_x[i])
            {
            outx_index = i;
            break;
            }
        }
//         gamma_x(outx_index) = gamma_x(len_gamma);
      gamma_x[outx_index] = gamma_x[len_gamma - 1];
//         gamma_x(len_gamma) = out_x;
      gamma_x[len_gamma - 1] = out_x[0]; // this doesn't do anything,
                                         // next line removes last
                                         // element of gamma_x
//         gamma_x = gamma_x(1:len_gamma-1);
      gamma_x.pop_back();
//         gamma_xk = gamma_x;
      gamma_xk = gamma_x;
//
//         rowi = outx_index; % ith row of A is swapped with last row (out_x)
      unsigned long rowi = outx_index;
//         colj = outx_index; % jth column of A is swapped with last column (out_lambda)
      unsigned long colj = outx_index;
//         AtgxAgx_ij = AtgxAgx;
      MatrixType AtgxAgx_ij = AtgxAgx;
//         temp_row = AtgxAgx_ij(rowi,:);
      MatrixType temp_row = AtgxAgx_ij.row(rowi);
//         AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
      AtgxAgx_ij.row(rowi) = AtgxAgx_ij.row(len_gamma - 1);
//         AtgxAgx_ij(len_gamma,:) = temp_row;
      AtgxAgx_ij.row(len_gamma - 1) = temp_row;
//         temp_col = AtgxAgx_ij(:,colj);
      MatrixType temp_col = AtgxAgx_ij.col(colj);
//         AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
      AtgxAgx_ij.col(colj) = AtgxAgx_ij.col(len_gamma - 1);
//         AtgxAgx_ij(:,len_gamma) = temp_col;
      AtgxAgx_ij.col(len_gamma - 1) = temp_col;
//         iAtgxAgx_ij = iAtgxAgx;
      MatrixType iAtgxAgx_ij = iAtgxAgx;
//         temp_row = iAtgxAgx_ij(colj,:);
      temp_row = iAtgxAgx_ij.row(colj);
//         iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
      iAtgxAgx_ij.row(colj) = iAtgxAgx_ij.row(len_gamma - 1);
//         iAtgxAgx_ij(len_gamma,:) = temp_row;
      iAtgxAgx_ij.row(len_gamma - 1) = temp_row;
//         temp_col = iAtgxAgx_ij(:,rowi);
      temp_col = iAtgxAgx_ij.col(rowi);
//         iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
      iAtgxAgx_ij.col(rowi) = iAtgxAgx_ij.col(len_gamma - 1);
//         iAtgxAgx_ij(:,len_gamma) = temp_col;
      iAtgxAgx_ij.col(len_gamma-1) = temp_col;
//
//         AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
      AtgxAgx = AtgxAgx_ij.block(0,0,len_gamma-1,len_gamma-1);
//         iAtgxAgx = update_inverse(AtgxAgx_ij, iAtgxAgx_ij,2);
      iAtgxAgx = AtgxAgx.inverse();
//         xk_1(out_x) = 0;
      for(unsigned int i = 0; i < out_x.size(); ++i)
        {
        xk_1(out_x[i],0) = 0.0;
        }
      }
//     else
    else
      {
//         % If an element is added to gamma_x
//         gamma_xk = [gamma_x; i_delta];
      gamma_xk = gamma_x;
      gamma_xk.push_back(i_delta);
//         new_x = i_delta;
      unsigned long new_x = i_delta;
//
//         AtgxAnx = A(:,gamma_x)'*A(:,new_x);
      MatrixType AtgxAnx;
      {
      MatrixType Atmp(A.rows(),gamma_x.size());
      for(unsigned int i = 0; i < gamma_x.size(); ++i)
        {
        Atmp.col(i) = A.col(gamma_x[i]);
        }
      AtgxAnx = Atmp.transpose() * A.col(new_x);
      }
//         AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];
      //
      // makes a new matrix of this form
      // | AtgxAgx  AtgxAnx |
      // | AtgxAnx' X       | where X = A(:,new_x)'*A(:,i_delta)
      // and this works great if AtgxAgx.rows() == AtgxAnx rows,
      // and AtgxAgx is square.
      MatrixType AtgxAgx_mod(AtgxAgx.rows() + 1, AtgxAgx.cols() + 1);
      {
      unsigned rows = AtgxAgx.rows(), cols = AtgxAgx.cols();
      AtgxAgx_mod.block(0,0,rows,cols) = AtgxAgx;
      AtgxAgx_mod.block(0,cols,AtgxAnx.rows(),AtgxAnx.cols()) = AtgxAnx;
      AtgxAgx_mod.block(rows,0,AtgxAnx.cols(),AtgxAnx.rows()) = AtgxAnx.transpose();
      AtgxAgx_mod(rows,cols) = A.col(new_x).transpose() * A.col(i_delta);
      }
//
//         AtgxAgx = AtgxAgx_mod;
      AtgxAgx = AtgxAgx_mod;
//         iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
      iAtgxAgx = AtgxAgx.inverse();
//         xk_1(i_delta) = 0;
      xk_1(i_delta) = 0;
//     end
      }
//
//     z_xk = zeros(N,1);
    z_xk = MatrixType::Zero(N,1);
//     z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
    for(unsigned int i = 0; i < gamma_xk.size(); ++i)
      {
      z_xk(gamma_xk[i],0) = -sign(Primal_constrk(gamma_xk[i],0));
      }
//     Primal_constrk([gamma_x]) = sign(Primal_constrk([gamma_x]))*epsilon;
    for(unsigned int i = 0; i < gamma_x.size(); ++i)
      {
      Primal_constrk(gamma_x[i],0) = sign(Primal_constrk(gamma_x[i])) * epsilon;
      }
// end
    }
// total_iter = iter;
// x_out = xk_1;
  MatrixType x_out = xk_1.transpose();
  return x_out;
}

// % BPDN_homotopy_function.m
// %
// % Solves the following basis pursuit denoising (BPDN) problem
// % min_x  \tau ||x||_1 + 1/2*||y-Ax||_2^2
// %
// % Inputs:
// % A - m x n measurement matrix
// % y - measurement vector
// % tau - final value of regularization parameter
// % maxiter - maximum number of homotopy iterations
// %
// % Outputs:
// % x_out - output for BPDN
// % gamma_x - support of the solution
// % total_iter - number of homotopy iterations taken by the solver
// % total_time - time taken by the solver
// %
// % Written by: Salman Asif, Georgia Tech
// % Email: sasif@ece.gatech.edu
// %
// %-------------------------------------------+
// % Copyright (c) 2007.  Muhammad Salman Asif
// %-------------------------------------------+
//
// function [x_out, gamma_x, total_iter, total_time] = BPDN_homotopy_function(A, y, tau, maxiter)
//
// t0 = cputime;
//
// N = size(A,2);
// K = size(A,1);
//
// % Initialization of primal and dual sign and support
// z_x = zeros(N,1);
// gamma_x = [];       % Primal support
//
// % Initial step
// Primal_constrk = -A'*y;
// [c i] = max(abs(Primal_constrk));
//
// gamma_xk = i;
//
// epsilon = c;
// xk_1 = zeros(N,1);
//
// z_x(gamma_xk) = -sign(Primal_constrk(gamma_xk));
// Primal_constrk(gamma_xk) = sign(Primal_constrk(gamma_xk))*epsilon;
//
// z_xk = z_x;
//
// % loop parameters
// done = 0;
// iter = 0;
// data_precision = eps;   % floating point precision
//
// old_delta = 0;
// out_x = [];
// count_delta_stop = 0;
//
// constraint_plots = 1;
//
// AtgxAgx = A(:,gamma_xk)'*A(:,gamma_xk);
// iAtgxAgx = inv(A(:,gamma_xk)'*A(:,gamma_xk));
//
// while iter < maxiter
//     iter = iter+1;
//     % warning('off','MATLAB:divideByZero')
//
//     gamma_x = gamma_xk;
//     z_x = z_xk;
//     x_k = xk_1;
//
//     %%%%%%%%%%%%%%%%%%%%%
//     %%%% update on x %%%%
//     %%%%%%%%%%%%%%%%%%%%%
//
//     % Update direction
//     %del_x = inv(A(:,gamma_x)'*A(:,gamma_x))*z_x(gamma_x);
//     del_x = iAtgxAgx*z_x(gamma_x);
//     del_x_vec = zeros(N,1);
//     del_x_vec(gamma_x) = del_x;
//
//     pk = Primal_constrk;
//     %dk = A'*(A*del_x_vec);
//     Agdelx = A(:,gamma_x)*del_x;
//     dk = A'*Agdelx;
//
//     %%% CONTROL THE MACHINE PRECISION ERROR AT EVERY OPERATION: LIKE BELOW.
//     pk_temp = Primal_constrk;
//     gammaL_temp = find(abs(abs(Primal_constrk)-epsilon)<min(epsilon,2*eps));
//     pk_temp(gammaL_temp) = sign(Primal_constrk(gammaL_temp))*epsilon;
//
//     xk_temp = x_k;
//     gammaX_temp = find(abs(x_k)<2*eps);
//     xk_temp(gammaX_temp) = 0;
//     %%%---
//
//     % Compute the step size
//     [i_delta, out_x, delta, chk_x] = update_primal(gamma_x, gamma_x, z_x,  xk_temp, del_x_vec, pk_temp, dk, epsilon, out_x);
//
//     if old_delta < 4*eps && delta < 4*eps
//         count_delta_stop = count_delta_stop + 1;
//     else
//         count_delta_stop = 0;
//     end
//     if count_delta_stop >= 500
//         %disp('stuck in some corner');
//         break;
//     end
//     old_delta = delta;
//
//     xk_1 = x_k+delta*del_x_vec;
//     Primal_constrk = pk+delta*dk;
//     epsilon_old = epsilon;
//     epsilon = epsilon-delta;
//
//     if epsilon <= tau;
//         xk_1 = x_k + (epsilon_old-tau)*del_x_vec;
//         total_time= cputime-t0;
//         break;
//     end
//
//     if chk_x == 1
//         % If an element is removed from gamma_x
//         gx_old = gamma_x;
//         len_gamma = length(gamma_x);
//
//         outx_index = find(gamma_x==out_x);
//         gamma_x(outx_index) = gamma_x(len_gamma);
//         gamma_x(len_gamma) = out_x;
//         gamma_x = gamma_x(1:len_gamma-1);
//         gamma_xk = gamma_x;
//
//         rowi = outx_index; % ith row of A is swapped with last row (out_x)
//         colj = outx_index; % jth column of A is swapped with last column (out_lambda)
//         AtgxAgx_ij = AtgxAgx;
//         temp_row = AtgxAgx_ij(rowi,:);
//         AtgxAgx_ij(rowi,:) = AtgxAgx_ij(len_gamma,:);
//         AtgxAgx_ij(len_gamma,:) = temp_row;
//         temp_col = AtgxAgx_ij(:,colj);
//         AtgxAgx_ij(:,colj) = AtgxAgx_ij(:,len_gamma);
//         AtgxAgx_ij(:,len_gamma) = temp_col;
//         iAtgxAgx_ij = iAtgxAgx;
//         temp_row = iAtgxAgx_ij(colj,:);
//         iAtgxAgx_ij(colj,:) = iAtgxAgx_ij(len_gamma,:);
//         iAtgxAgx_ij(len_gamma,:) = temp_row;
//         temp_col = iAtgxAgx_ij(:,rowi);
//         iAtgxAgx_ij(:,rowi) = iAtgxAgx_ij(:,len_gamma);
//         iAtgxAgx_ij(:,len_gamma) = temp_col;
//
//         AtgxAgx = AtgxAgx_ij(1:len_gamma-1,1:len_gamma-1);
//         iAtgxAgx = update_inverse(AtgxAgx_ij, iAtgxAgx_ij,2);
//         xk_1(out_x) = 0;
//     else
//         % If an element is added to gamma_x
//         gamma_xk = [gamma_x; i_delta];
//         new_x = i_delta;
//
//         AtgxAnx = A(:,gamma_x)'*A(:,new_x);
//         AtgxAgx_mod = [AtgxAgx AtgxAnx; AtgxAnx' A(:,new_x)'*A(:,i_delta)];
//
//         AtgxAgx = AtgxAgx_mod;
//         iAtgxAgx = update_inverse(AtgxAgx, iAtgxAgx,1);
//         xk_1(i_delta) = 0;
//     end
//
//     z_xk = zeros(N,1);
//     z_xk(gamma_xk) = -sign(Primal_constrk(gamma_xk));
//     Primal_constrk([gamma_x]) = sign(Primal_constrk([gamma_x]))*epsilon;
// end
// total_iter = iter;
// x_out = xk_1;
