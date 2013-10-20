// A common file to define types
// This will ease the transition
// using new types in eigen.
#ifndef __ukf_type_h__
#define __ukf_type_h__

#include "Eigen/Dense"

typedef double ukfPrecisionType;

/** Short hand for the state vector */
#include <vector>
#if 0
typedef std::vector<ukfPrecisionType> State;
#else
typedef std::vector<ukfPrecisionType> stdVecState;
typedef Eigen::VectorXd  State;
#endif

typedef Eigen::VectorXd ukfVectorType;
typedef Eigen::MatrixXd ukfMatrixType;

/**
 * \struct vec_t
 * \brief Defines simple 3D vector and matrices. 3
*/
typedef Eigen::Vector3d  vec_t;
typedef std::vector<vec_t,Eigen::aligned_allocator<vec_t> > stdVec_t;
typedef std::vector<ukfVectorType,Eigen::aligned_allocator<ukfVectorType> > stdEigVec_t;

/**
 * \struct mat_t
 * \brief Defines a 3-by-3 Matrix
*/
typedef Eigen::Matrix<ukfPrecisionType,3,3> mat_t;
typedef std::vector<mat_t,Eigen::aligned_allocator<mat_t> > stdMat_t;

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
