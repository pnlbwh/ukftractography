#ifndef __CompressedSensing_h
#define __CompressedSensing_h
#include "ukf_types.h"
#include "itkImage.h"
typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
typedef std::vector<MatrixType> MatrixVector;
typedef Eigen::Matrix<ukfPrecisionType,Eigen::Dynamic,1> VectorType;
typedef itk::Image<unsigned char,3> MaskImageType;
#endif // __CompressedSensing_h
