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

/**
 * \struct vec_t
 * \brief Defines simple 3D vector and matrices. 3
*/
typedef Eigen::Vector3d  vec_t;

/**
 * \struct mat_t
 * \brief Defines a 3-by-3 Matrix
*/
typedef Eigen::Matrix<double,3,3> mat_t;

typedef vnl_vector<double> ukfVectorType;
typedef vnl_matrix<double> ukfMatrixType;

#endif // __ukf_type_h__
