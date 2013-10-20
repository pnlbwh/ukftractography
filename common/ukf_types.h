// A common file to define types
// This will ease the transition
// using new types in eigen.
#ifndef __ukf_type_h__
#define __ukf_type_h__

#include "Eigen/Dense"

typedef double ukfPrecisionType;
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
