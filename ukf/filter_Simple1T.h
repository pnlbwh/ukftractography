#ifndef SIMPLE1T_H__
#define SIMPLE1T_H__
#include "filter_model.h"
/**
 * \struct Simple1T
 * \brief Simple 1-Tensor model
 *
 * Model describing 1-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal)
*/
class Simple1T : public FilterModel
{
public:
  Simple1T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
           bool constrained)
    : FilterModel(5, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(3, 3) = _Q(4, 4) = ql;
  }

  virtual ~Simple1T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor1T(const State& x, vec3_t& m, vec3_t& l);

  const ukfPrecisionType _lambda_min;
};

#endif // SIMPLE1T_H__
