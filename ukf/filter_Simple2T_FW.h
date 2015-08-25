#ifndef SIMPLE2T_FW__H_
#define SIMPLE2T_FW__H_
#include "filter_model.h"

/**
 * \struct Simple2T_FW
 * \brief Simple 2-Tensor model with free water
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
 * and free water estimation
*/
class Simple2T_FW : public FilterModel
{
public:
  Simple2T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs,
              const ukfVectorType& weights_on_tensors,
              bool constrained, const ukfPrecisionType diff_fw)
    : FilterModel(11, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw) )
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
    _Q(10, 10) = qw; // noise for weights

    _D.resize(11, 6);
    _D.setConstant(ukfZero);

    _d.resize(6);

    // Setting the constraints according to D'*x >= -d
    _D(10, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(10, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2)  = 1;
    _d(2) = 0;  // l11 >= 0
    _D(4, 3)  = 1;
    _d(3) = 0;  // l21 >= 0
    _D(8, 4)  = 1;
    _d(4) = 0;  // l12 >= 0
    _D(9, 5)  = 1;
    _d(5) = 0;  // l22 >= 0
  }

  virtual ~Simple2T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;
};
#endif// SIMPLE2T_FW__H_
