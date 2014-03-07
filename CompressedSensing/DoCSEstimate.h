#ifndef __DoCSEstimate_h
#define __DoCSEstimate_h
#include "CompressedSensing.h"
#include "NrrdFile.h"
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
  void step2(DWIVectorImageType::Pointer &inputImage,double myu, double tol,DWIVectorImageType::Pointer &rval);
  MatrixType step1(DWIVectorImageType::Pointer &inputImage,
                   MatrixType &A,
                   double lmd,
                   unsigned NIT,
                   std::vector<unsigned char> &id);
  void tvdenoise3(DWIVectorImageType::Pointer &inputImage,
                  unsigned int gradientIndex,double lambda,
                  double tol,
                  DWIVectorImageType::Pointer &target);
  void ToVecImage(const B0AvgImageType *fromImage, unsigned int gradientIndex, DWIVectorImageType::Pointer &target);
  B0AvgImageType *FromVecImage(DWIVectorImageType::Pointer inputImage,unsigned gradientIndex);

  // MatrixType BPDN_HOMOTOPY_function(MatrixType &A,MatrixType &SqueezeS,double tau, unsigned int maxiter);
  // void update_primal(std::vector<unsigned long> & gamma_x, // current support of x
  //                    std::vector<unsigned long> & gamma_lambda, // current support of
  //                    // lambda
  //                    MatrixType &                 z_x, // sign sequence of x
  //                    MatrixType &                 x_k, // sign sequence of lambda
  //                    MatrixType &                 del_x_vec, // primal update direction
  //                    MatrixType &                 pk,       //
  //                    MatrixType &                 dk,       //
  //                    double                       epsilon, // current value of epsilon
  //                    std::vector<unsigned long> & out_lambda, // element removed from
  //                                                             // support of lambda in
  //                                                             // previous step if any
  //                                                             // OUTPUTS
  //                    unsigned long &              i_delta, // index corresponding to
  //                                                          // newly active primal
  //                                                          // constraint (new_lambda)
  //                    std::vector<unsigned long> & out_x, // element in x shrunk to zero;
  //                    double  &                    delta, // primal step size
  //                    unsigned &                   chk_x); // 1 an element is removed
  //                                                        // from support of x
  //                                                        // 0 a new element enters
  //                                                        // the support of lambd

  void Reshape(const DWIVectorImageType::Pointer &src,
               unsigned rows, unsigned columns,
               MatrixType &outMatrix);
  void Reshape(const MatrixType &src,
               DWIVectorImageType::Pointer templateImage,
               DWIVectorImageType::Pointer &outMatrix);
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
