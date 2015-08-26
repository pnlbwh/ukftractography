#ifndef NODDI1F_H__
#define NODDI1F_H__
#include "filter_model.h"

/**
 * \struct NODDI1F
 * \brief NOODI 1 fiber model
 *
 * Model describing 1-fiber tractography with the NODDI representation
*/
class NODDI1F : public FilterModel
{
public:
  NODDI1F(ukfPrecisionType qs, ukfPrecisionType qkappa, ukfPrecisionType qvic, ukfPrecisionType rs,
          const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(6, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = qs;
    _Q(4, 4) = qkappa;
    _Q(3, 3) = _Q(5, 5) = qvic;

    _D.resize(6, 5);
    _D.setConstant(ukfZero);

    _d.resize(5);

    // Setting the constraints according to D'*x >= -d
    _D(5, 0) = -1;
    _d(0) = 1;  // Viso <= 1
    _D(5, 1) = 1;
    _d(1) = 0;  // Viso >= 0

    _D(3, 2) = 1;
    _d(2) = 0;  // Vic >= 0
    _D(3, 3) = -1;
    _d(3) = 1;  // Vic <= 1

    _D(4, 4) = 1;
    _d(4) = 0;  // kappa >= 0

  }

  virtual ~NODDI1F()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  const ukfPrecisionType _lambda_min;
};
#endif // NODDI1F_H__
