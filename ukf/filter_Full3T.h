#ifndef FULL3T_H__
#define FULL3T_H__
#include "filter_model.h"

/**
 * \struct Full3T
 * \brief Full 3-Tensor model with free water
 *
 * Model describing 3-tensor tractography with the full tensor representation
 * (3 angles, 3 eigenvalues).
*/
class Full3T : public FilterModel
{
public:
  Full3T(ukfPrecisionType qs, ukfPrecisionType ql, ukfPrecisionType rs, const ukfVectorType& weights_on_tensors,
         bool constrained)
    : FilterModel(18, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(6, 6) = _Q(7, 7) = _Q(8, 8) = qs;
    _Q(12, 12) = _Q(13, 13) = _Q(14, 14) = qs;

    _Q(3, 3) = _Q(4, 4) = _Q(5, 5) = ql;
    _Q(9, 9) = _Q(10, 10) = _Q(11, 11) = ql;
    _Q(15, 15) = _Q(16, 16) = _Q(17, 17) = ql;
  }

  virtual ~Full3T()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  virtual void State2Tensor3T(const State& x, const vec3_t& old_m, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2,
                              vec3_t& m3, vec3_t& l3);

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

};
#endif // FULL3T_H__
