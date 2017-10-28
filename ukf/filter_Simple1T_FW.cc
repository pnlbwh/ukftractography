#include "filter_Simple1T_FW.h"

// Functions for 1-tensor simple model.
void Simple1T_FW::F(ukfMatrixType& X) const
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

    // Clamp lambdas.
//     X(3, i) = std::max(X(3, i), ukfZero); // 0 because for free water lambdas
// are constrained to be > 0
//     X(4, i) = std::max(X(4, i), ukfZero); // because imprecission of qp they
// are sometimes slightly negative on the order
// of e-16
    X(3, i) = CheckZero(X(3, i) );
    X(4, i) = CheckZero(X(4, i) );

    X(5, i) = CheckZero(X(5, i) ); // fw

    // DEBUGGING
    if( X(3, i) < 0 || X(4, i) < 0 )
      {
      std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
      }
    }
}

void Simple1T_FW::H(const ukfMatrixType& X,
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

  const stdVec_t&       gradients = _signal_data->gradients();
  const ukfVectorType & b     = _signal_data->GetBValues();
  diagmat3_t            lambdas;
  for( unsigned int i = 0; i < X.cols(); ++i )
    {
    // Normalize direction.
    vec3_t m;
    initNormalized(m, X(0, i), X(1, i), X(2, i) );

    // Clamp lambdas.
    lambdas.diagonal()[0] = CheckZero(X(3, i) );
    lambdas.diagonal()[1] = CheckZero(X(4, i) );
    lambdas.diagonal()[2] = lambdas.diagonal()[1];
    if( lambdas.diagonal()[0] < 0 ||  lambdas.diagonal()[1]  < 0 )
      {
      std::cout << "Warning: eigenvalues became negative l1= " << lambdas.diagonal()[0] << "l2= "
                << lambdas.diagonal()[1] << std::endl;
      }

    // get weight from state
    const ukfPrecisionType w = CheckZero( X(5, i) );
    // FOR DEBUGGING :
    if( w < 0 - 1.0e-5 )
      {
      std::cout << "Negative weight! w= " << w << "\n";
      std::cout << X << "\ni: " << i;
      throw;
      }
    if( w > 1 + 1.0e-5 )
      {
      std::cout << "Weight > 1 => negative free water! w= " << w << "\n";
      std::cout << X << "\ni: " << i;
      throw;
      }

    // Flip if necessary.
    // Why is that???
    if( m[0] < 0 )
      {
      m = -m;
      }

    // Calculate diffusion matrix.
    const mat33_t & local_D = diffusion(m, lambdas); // l3 == l2
    // Reconstruct signal.
    for( int j = 0; j < _signal_dim; ++j )
      {
      const vec3_t& u = gradients[j];
      Y(j, i) =  (w) * std::exp(-b[j] * u.dot(local_D * u) )
        + (1 - w) * std::exp(-b[j] * u.dot(m_D_iso * u) );
      }
    }
}

void Simple1T_FW::State2Tensor1T(const State& x, vec3_t& m, vec3_t& l)
{
  // Orientation.
  initNormalized(m, x[0], x[1], x[2]);

  l[1] = std::max(x[4], ukfZero );

  l[0] = CheckZero(x[3]);
  l[1] = CheckZero(x[4]);
  if( l[0] < 0 || l[1] < 0 )
    {
    std::cout << "Warning: eigenvalues became negative " << __LINE__ << " " << __FILE__ << std::endl;
    }
  l[2] = l[1];
}
