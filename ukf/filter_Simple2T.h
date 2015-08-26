#ifndef SIMPLE2T_H__
#define SIMPLE2T_H__
#include "filter_model.h"
/**
 * \struct Simple2T
 * \brief Simple 2-Tensor model
 *
 * Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal).
*/
class Simple2T : public FilterModel
{
public:
  Simple2T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
           bool constrained)
    : FilterModel(10, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(8, 8) = _Q(9, 9) = ql;
  }

  virtual ~Simple2T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor2T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

};

#endif // SIMPLE2T_H__
