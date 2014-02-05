#ifndef __DoCSEstimate_h
#define __DoCSEstimate_h
#include "CompressedSensing.h"
#include "nrrdIO.h"
#include <vector>
#include <cmath>
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "icosahedron3.h"
#include "SphericalHarmonics.h"

class DoCSEstimate
{
public:
  typedef itk::Image<double,3> B0AvgImageType;
  typedef NrrdFile::ImageType  DWIVectorImageType;

  DoCSEstimate(NrrdFile &nrrdFile,MaskImageType::Pointer &maskImage,MatrixType &newGradients);
  bool AverageB0AndExtractIntensity();
  bool Compute();
  void OptimBand(double B,const MatrixType &D,const MatrixType &ico3,double  &rho, double  &p);
private:
  NrrdFile &                    m_NrrdFile;
  MaskImageType::Pointer        m_MaskImage;
  DWIVectorImageType::PointType m_SpaceOrigin;
  MatrixType                    m_NewGradients;
  double                        m_BValue;
  MatrixType                    m_SpaceDirections;
  MatrixType                    m_VoxelLatticeAlignedGradientDirections;
  DWIVectorImageType::Pointer   m_IntensityData;
  B0AvgImageType::Pointer       m_AverageB0;
  const double                  m_Pi;
};
#endif // __DoCSEstimate_h
