#ifndef SIMPLE2BIEXP_FW_H__
#define SIMPLE2BIEXP_FW_H__
#include "filter_model.h"



/**
 * \struct Simple2T_BiExp_FW
 * \brief Simple 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
class Simple2T_BiExp_FW : public FilterModel
  {
public:
  Simple2T_BiExp_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
              bool constrained, const ukfPrecisionType diff_fw, ukfPrecisionType qt)
    : FilterModel(15, rs, weights_on_tensors, constrained), _lambda_min_fast_diffusion(1.0), _lambda_min_slow_diffusion(0.1), _lambda_max_diffusion(3000), _w_fast_diffusion(0.7), m_D_iso(SetIdentityScaled(diff_fw))
  {

    // 15 = size(X, 'column')
    // X = [x10, x11, x12, l11, l12, l13, l14, x20, x21, x22, l21, l22, l23, l24, w]' 
    //     [0  , 1  , 2  , 3  , 4  , 5  , 6  , 7  , 8  , 9  , 10 , 11 , 12 , 13 , 14]
    // There are two tensors in this model, and we have a simple bi-exponential model, which means 2 lambdas per exponential (simple model), and 2 exponentials per tensor (bi-exponential)
    // 1-w is the free water weight


    // Q is the injected process noise bias. It is a diagonal matrix. It is used in the state transfer function.

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(7, 7) = _Q(8, 8) = _Q(9, 9) = qs;                                   // noise for the main direction
    _Q(3, 3) = _Q(4, 4) =  _Q(10, 10) = _Q(11, 11) = ql;                                                    // noise for the lambdas (fast diffusion)
    _Q(5, 5) = _Q(6, 6) =  _Q(12, 12) = _Q(13, 13) = qt;                                                    // noise for the lambdas (slow diffusion)
    _Q(14, 14) = qw;                                                                                        // noise for the free water weight

    // D is the constraint matrix. 
    // The 1st dim of D is used to select which component of x we want to constraint, the 2nd is for the inequality (*-1 -> reverse the inequality)
    // d is the contraint value
    // D'*x >= -d
    
    // 18 constraints for the 15 dimensions of the state
    _D.resize(15, 18);
    _D.setConstant(ukfZero);

    _d.resize(18);

    // Setting the constraints according to D'*x >= -d
    
    // Free water
    _D(14, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(14, 1) = 1;
    _d(1) = 0;  // w >= 0

    
    // Tensor 1 (minimal values)
    _D(3, 2)  = 1;
    _d(2) = _lambda_min_fast_diffusion;  // l11 >= min
    _D(4, 3)  = 1;
    _d(3) = _lambda_min_fast_diffusion;  // l12 >= min
    _D(5, 4)  = 1;
    _d(4) = _lambda_min_slow_diffusion;  // l13 >= min
    _D(6, 5)  = 1;
    _d(5) = _lambda_min_slow_diffusion;  // l14 >= min
    
    
    // Tensor 2 (minimal values)
    _D(10, 6)  = 1;
    _d(6) = _lambda_min_fast_diffusion;  // l21 >= min
    _D(11, 7)  = 1;
    _d(7) = _lambda_min_fast_diffusion;  // l22 >= min
    _D(12, 8)  = 1;
    _d(8) = _lambda_min_slow_diffusion;  // l23 >= min
    _D(13, 9)  = 1;
    _d(9) = _lambda_min_slow_diffusion;  // l21 >= min
    
    // Tensor 1 (maximal values)
    _D(3, 10) = -1;                      // l11 <= max
    _d(10) = _lambda_max_diffusion;
    _D(4, 11) = -1;                      // l12 <= max
    _d(11) = _lambda_max_diffusion;
    _D(5, 12) = -1;                      // l13 <= max
    _d(12) = _lambda_max_diffusion;
    _D(6, 13) = -1;                      // l14 <= max
    _d(13) = _lambda_max_diffusion;
    
    // Tensor 2 (maximal values)
    _D(10, 14) = -1;                      // l21 <= max
    _d(14) = _lambda_max_diffusion;
    _D(11, 15) = -1;                      // l22 <= max
    _d(15) = _lambda_max_diffusion;
    _D(12, 16) = -1;                      // l23 <= max
    _d(16) = _lambda_max_diffusion;
    _D(13, 17) = -1;                      // l24 <= max
    _d(17) = _lambda_max_diffusion;
    

  }

  virtual ~Simple2T_BiExp_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l11, vec3_t& m2, vec3_t& l21);

  /** The minimum/maximum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min_fast_diffusion;
  const ukfPrecisionType _lambda_min_slow_diffusion;
  const ukfPrecisionType _lambda_max_diffusion;
  
  /** The weight fraction of the fast diffusion component */
  const ukfPrecisionType _w_fast_diffusion;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;
  };


#endif  // SIMPLE2BIEXP_FW_H__