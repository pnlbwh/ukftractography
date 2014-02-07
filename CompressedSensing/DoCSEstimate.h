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
  typedef itk::Image<double,3>         B0AvgImageType;
  typedef NrrdFile::ImageType          DWIVectorImageType;
  typedef DWIVectorImageType::SizeType ImageSizeType;
  typedef ImageSizeType::SizeValueType ImageSizeValueType;

  DoCSEstimate(NrrdFile &nrrdFile,MaskImageType::Pointer &maskImage,MatrixType &newGradients);
  bool AverageB0AndExtractIntensity();
  bool Compute();
  void OptimBand(double B,const MatrixType &D,const MatrixType &ico3,double  &rho, double  &p);
  DWIVectorImageType *step2(double myu, double tol);
  MatrixType step1(MatrixType &A,double lmd,unsigned NIT,std::vector<unsigned int> &id);
  void tvdenoise3(unsigned int gradientIndex,double lambda,double tol,DWIVectorImageType::Pointer &target);
  void ToVecImage(const B0AvgImageType *fromImage, unsigned int gradientIndex, DWIVectorImageType::Pointer &target);
  B0AvgImageType *FromVecImage(unsigned gradientIndex);
  MatrixType BPDN_HOMOTOPY_function(MatrixType &A,MatrixType &SqueezeS,double lmd, unsigned int NIT);
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
