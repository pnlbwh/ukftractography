// A common file to define types
// This will ease the transition
// using new types in eigen.
#ifndef __ukf_type_h__
#define __ukf_type_h__

#include "Eigen/Dense"
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

/** Short hand for the state vector */
#if 1
#include <vector>
typedef std::vector<double> State;
#else
typedef Eigen::VectorXd  State;
#endif

#if 0
typedef vnl_vector<double> ukfVectorType;
typedef vnl_matrix<double> ukfMatrixType;
#else
typedef Eigen::VectorXd ukfVectorType;
typedef Eigen::MatrixXd ukfMatrixType;
#endif

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
typedef Eigen::Matrix<double,3,3> mat_t;
typedef std::vector<mat_t,Eigen::aligned_allocator<mat_t> > stdMat_t;

// HACK:  A quick conversion code for testing conversions between libraries.
inline Eigen::MatrixXd ToMatrixXd( const vnl_matrix<double> & in)
{
  Eigen::MatrixXd out;
  out.resize(in.rows(),in.cols());
  for(int r=0; r< static_cast<int>(in.rows()); ++r)
    {
    for(int c=0; c< static_cast<int>(in.cols()); ++c)
      {
      out(r,c)=in[r][c];
      }
    }
  return out;
}
inline vnl_matrix<double> ToVNL( const Eigen::MatrixXd & in)
{
  vnl_matrix<double> out;
  out.set_size(in.rows(),in.cols());
  for(int r=0; r< in.rows(); ++r)
    {
    for(int c=0; c< in.cols(); ++c)
      {
      out[r][c]=in(r,c);
      }
    }
  return out;
}

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
