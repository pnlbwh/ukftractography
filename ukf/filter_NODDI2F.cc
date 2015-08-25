#include "filter_NODDI2F.h"
#include "utilities.h"

// Functions for 2-fiber noddi model.
void NODDI2F::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  // Clamp lambdas.
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vectors.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;

    norm_inv = ukfZero; // 1e-16;
    norm_inv += X(5, i) * X(5, i);
    norm_inv += X(6, i) * X(6, i);
    norm_inv += X(7, i) * X(7, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(5, i) *= norm_inv;
    X(6, i) *= norm_inv;
    X(7, i) *= norm_inv;
    }
}

void NODDI2F::H(const  ukfMatrixType& X,
                ukfMatrixType& Y) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         (X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          X.cols() == 1) );
  assert(Y.rows() == static_cast<unsigned int>(_signal_dim) &&
         (Y.cols() == static_cast<unsigned int>(2 * _state_dim + 1) ||
          Y.cols() == 1) );
  assert(_signal_data);

  ukfPrecisionType dPar, dIso;
  dPar = 0.0000000017;
  dIso = 0.000000003;
  const stdVec_t& gradients = _signal_data->gradients();
  ukfVectorType   gradientStrength, pulseSeparation;
  createProtocol(_signal_data->GetBValues(), gradientStrength, pulseSeparation);
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize directions.
    vec3_t m1;
    initNormalized(m1, X(0, i), X(1, i), X(2, i) );
    vec3_t m2;
    initNormalized(m2, X(5, i), X(6, i), X(7, i) );

    // Flip if necessary.
    if( m1[0] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] < 0 )
      {
      m2 = -m2;
      }

    ukfVectorType          Eec1, Eic1, Eec2, Eic2, Eiso;
    const ukfPrecisionType Vic1 = X(3, i);
    const ukfPrecisionType kappa1 = X(4, i);
    const ukfPrecisionType Vic2 = X(8, i);
    const ukfPrecisionType kappa2 = X(9, i);
    const ukfPrecisionType Viso = X(10, i);

    ExtraCelluarModel(dPar, Vic1, kappa1, gradientStrength,
                      pulseSeparation, gradients, m1, Eec1);
    IntraCelluarModel(dPar, kappa1, gradientStrength, pulseSeparation,
                      gradients, m1, Eic1 );

    ExtraCelluarModel(dPar, Vic2, kappa2, gradientStrength,
                      pulseSeparation, gradients, m2, Eec2);
    IntraCelluarModel(dPar, kappa2, gradientStrength, pulseSeparation,
                      gradients, m2, Eic2 );
    IsoModel(dIso, gradientStrength, pulseSeparation, Eiso);
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      Y(j, i) = Viso * Eiso[j] + (1 - Viso) * ( (Vic1 * Eic1[j] + (1 - Vic1) * Eec1[j]) * weights_on_tensors_[0]
                                                + (Vic2 * Eic2[j] + (1 - Vic2) * Eec2[j]) * weights_on_tensors_[1]);
      }
    }
}
