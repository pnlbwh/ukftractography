#ifndef FULL1T_H__
#define FULL1T_H__
#include "filter_model.h"

/**
 * \struct Full1T
 * \brief Full 1-Tensor model without free water
 *
 * Model describing 1-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
class Full1T : public FilterModel
{
public:
  Full1T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
         bool constrained)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
  }

  virtual ~Full1T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;
};

#endif // FULL1T_H__
