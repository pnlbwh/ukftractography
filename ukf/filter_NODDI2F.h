#ifndef NODDI2F_H__
#define NODDI2F_H__
#include "filter_model.h"
/**
* \struct Simple2T
* \brief Simple 2-Tensor model
*
* Model describing 2-tensor tractography with the simplified tensor representation (two minor eigenvalues are equal).
*/
class NODDI2F : public FilterModel
{
public:
  NODDI2F(ukfPrecisionType qs, ukfPrecisionType qkappa, ukfPrecisionType qvic, ukfPrecisionType rs,
          const ukfVectorType& weights_on_tensors, bool constrained)
    : FilterModel(11, rs, weights_on_tensors, constrained), _lambda_min(100.0)
  {
    _Q(0, 0) = _Q(1, 1) = _Q(2, 2) = _Q(5, 5) = _Q(6, 6) = _Q(7, 7) = qs;
    _Q(4, 4) = _Q(9, 9) = qkappa;
    _Q(3, 3) = _Q(8, 8) = _Q(10, 10) = qvic; // noise for weights

    _D.resize(11, 8);
    _D.setConstant(ukfZero);

    _d.resize(8);

    // Setting the constraints according to D'*x >= -d
    _D(10, 0) = -1;
    _d(0) = 1;  // Viso <= 1
    _D(10, 1) = 1;
    _d(1) = 0;  // Viso >= 0

    _D(3, 2) = 1;
    _d(2) = 0;  // Vic >= 0
    _D(3, 3) = -1;
    _d(3) = 1;  // Vic <= 1

    _D(4, 4) = 1;
    _d(4) = 0;  // kappa >= 0

    _D(8, 5) = 1;
    _d(5) = 0;  // Vic >= 0
    _D(8, 6) = -1;
    _d(6) = 1;  // Vic <= 1

    _D(9, 7) = 1;
    _d(7) = 0;  // kappa >= 0

  }

  virtual ~NODDI2F()
  {
  }

  virtual void F(ukfMatrixType& X) const;

  virtual void H(const  ukfMatrixType& X, ukfMatrixType& Y) const;

  /** The minimum value of the eigenvalues. Clamped in each step */
  const ukfPrecisionType _lambda_min;

};
#endif // NODDI2F_H__
