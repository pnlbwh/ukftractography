#include "filter_NODDI1F.h"
#include "utilities.h"

// Functions for 1-fiber NODDI model.
void NODDI1F::F(ukfMatrixType& X) const
{
  assert(_signal_dim > 0);
  assert(X.rows() == static_cast<unsigned int>(_state_dim) &&
         X.cols() == static_cast<unsigned int>(2 * _state_dim + 1) );
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize the direction vector.
    ukfPrecisionType norm_inv = ukfZero; // 1e-16;
    norm_inv += X(0, i) * X(0, i);
    norm_inv += X(1, i) * X(1, i);
    norm_inv += X(2, i) * X(2, i);

    norm_inv = ukfOne / sqrt(norm_inv);
    X(0, i) *= norm_inv;
    X(1, i) *= norm_inv;
    X(2, i) *= norm_inv;
    }
}

void NODDI1F::H(const  ukfMatrixType& X, ukfMatrixType& Y) const
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
    ukfVectorType Eec, Eic, Eiso;
    // Normalize direction.
    vec3_t m;
    initNormalized(m, X(0, i), X(1, i), X(2, i) );

    // Clamp lambdas.
    const ukfPrecisionType Vic = X(3, i);
    const ukfPrecisionType kappa = X(4, i);
    const ukfPrecisionType Viso = X(5, i);
    ExtraCelluarModel(dPar, Vic, kappa, gradientStrength,
                      pulseSeparation, gradients, m, Eec);
    IntraCelluarModel(dPar, kappa, gradientStrength, pulseSeparation,
                      gradients, m, Eic);
    IsoModel(dIso, gradientStrength, pulseSeparation, Eiso);
    assert(Eic.size() == _signal_dim);
    for( int j = 0; j < _signal_dim; ++j )
      {
      Y(j, i) = Viso * Eiso[j] + (1 - Viso) * (Vic * Eic[j] + (1 - Vic) * Eec[j]);
      }
    }
}
