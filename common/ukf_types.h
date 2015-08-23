// A common file to define types
// This will ease the transition
// using new types in eigen.
#ifndef __ukf_type_h__
#define __ukf_type_h__

#include "Eigen/Dense"

//#define UKF_USE_FLOAT 1
#if UKF_USE_FLOAT
typedef float ukfPrecisionType;
#else
typedef double ukfPrecisionType;
#endif

static const ukfPrecisionType GLOBAL_TENSOR_PACK_VALUE=1e6;
static const ukfPrecisionType GLOBAL_TENSOR_UNPACK_VALUE=1e-6;

static const ukfPrecisionType ukfZero(static_cast<ukfPrecisionType>(0.0));
static const ukfPrecisionType ukfOne(static_cast<ukfPrecisionType>(1.0));
static const ukfPrecisionType ukfHalf(static_cast<ukfPrecisionType>(0.5));

typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,1> ukfVectorType;
typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,Eigen::Dynamic> ukfMatrixType;
typedef Eigen::Matrix<ukfPrecisionType,3,1> vec3_t;
typedef Eigen::Matrix<ukfPrecisionType,3,3> mat33_t;


/** Short hand for the state vector */
#include <vector>
typedef std::vector<ukfPrecisionType> stdVecState;
typedef ukfVectorType  State;

/**
 * \struct vec3_t
 * \brief Defines simple 3D vector and matrices. 3
*/
typedef std::vector<vec3_t,Eigen::aligned_allocator<vec3_t> > stdVec_t;
typedef std::vector<ukfVectorType,Eigen::aligned_allocator<ukfVectorType> > stdEigVec_t;
typedef std::vector<mat33_t,Eigen::aligned_allocator<mat33_t> > stdMat_t;

template <typename TIn, typename TOut>
TOut ConvertVector(const TIn &in)
{
  TOut out(in.size());
  //std::copy(in.begin(),in.end(),out.begin());
  for (unsigned int i =0 ; i < in.size(); ++i )
  {
    out[i] = in[i];
  }
  return out;
}

#endif // __ukf_type_h__
