#ifndef FULL1T_FW_H__
#define FULL1T_FW_H__
#include "filter_model.h"

/**
 * \struct Full1T_FW
 * \brief Full 1-Tensor model with free water
 *
 * Model describing 1-tensor tractography with the full tensor representation (3 angles, 3 eigenvalues)
 * and free water estimation.
*/
class Full1T_FW : public FilterModel
{
public:
  Full1T_FW(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType qw, ukfPrecisionType rs,
            const ukfVectorType& weights_on_tensors, bool constrained,
            const ukfPrecisionType diff_fw)
    : FilterModel(7, rs, weights_on_tensors, constrained), _lambda_min(100.0), m_D_iso(SetIdentityScaled(diff_fw) )
  {
#if 0
    m_D_iso << diff_fw, 0, 0,
      0, diff_fw, 0,
      0, 0, diff_fw;
#endif

    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
    _Q(6, 6) = qw;

    _D.resize(7, 5);
    _D.setConstant(ukfZero);

    _d.resize(5);

    // Setting the constraints according to D'*x >= -d
    _D(6, 0) = -1;
    _d(0) = 1;  // w <= 1
    _D(6, 1) = 1;
    _d(1) = 0;  // w >= 0
    _D(3, 2) = 1;
    _d(2) = 0;  // l1 >= 0
    _D(4, 3) = 1;
    _d(3) = 0;  // l2 >= 0
    _D(5, 4) = 1;
    _d(4) = 0;  // l3 >= 0

  }

  virtual ~Full1T_FW()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

  /** apparent diffusion coefficient of free water */
  mat33_t m_D_iso;

};
#endif // FULL1T_FW_H__
