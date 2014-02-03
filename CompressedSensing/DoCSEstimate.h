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

  DoCSEstimate(NrrdFile &nrrdFile,MaskImageType::Pointer &maskImage,MatrixType &newGradients) :
    m_NrrdFile(nrrdFile), m_MaskImage(maskImage), m_Gradients(newGradients)
    {
      this->m_SpaceOrigin = nrrdFile.GetImage()->GetOrigin();
      this->m_B0 = nrrdFile->GetBValue();
      DWIVectorImageType::DirectionType dir = nrrdFile.GetImage()->GetDirection();
      this->m_SpaceDirections.conservativeResize(3,3);
      for(unsigned int i = 0; i < 3; ++i)
        {
        for(unsigned int j = 0; j < 3; ++j)
          {
          this->m_SpaceDirections(i,j) = dir[i][j];
          }
        }
    }
  bool AverageB0AndExtractIntensity()
    {
      const double eps(2.2204e-16);

      MatrixType gradients = this->m_NrrdFile.GetGradients();
      std::vector<unsigned> b0Indices;
      std::vector<unsigned> gradientIndices;

      for(unsigned i = 0; i < gradients.rows(); ++i)
        {
        double norm = gradients.row(i).norm();
        if(norm < eps)
          {
          b0Indices.push_back(i);
          }
        else
          {
          gradientIndices.push_back(i);
          }
        }

      const unsigned int numB0(b0Indices.size());
      if(numB0 == 0)
        {
        std::cerr << "input data doesn't have a b0 image" << std::endl;
        return false;
        }

      // allocate b0 average
      DWIVectorImageType::Pointer originalNrrd = this->m_NrrdFile.GetImage();
      DWIVectorImageType::RegionType region = originalNrrd->GetLargestPossibleRegion();
      this->m_AverageB0 = B0AvgImageType::New();
      this->m_AverageB0->SetRegions(region);
      this->m_AverageB0->Allocate();
      this->m_AverageB0->FillBuffer(0.0);

      itk::ImageRegionConstIterator<DWIVectorImageType> b0It(originalNrrd,region);
      itk::ImageRegionIterator<B0AvgImageType> avgIt(this->m_AverageB0,region);

      for(b0It.Begin(), avgIt.Begin(); !b0It.IsAtEnd() && !avgIt.IsAtEnd(); ++b0It, ++avgIt)
        {
        const DWIVectorImageType::PixelType val = b0It.Get();
        B0AvgImageType::PixelType avg(0.0);
        for(unsigned int i = 0; i < numB0; ++i)
          {
          avg += val[b0Indices[i]];
          }
        avg /= static_cast<double>(numB0);
        if(std::fabs(avg) <= eps)
          {
          avg = 1.0;
          }
        avgIt.Set(avg);
        }
      //
      // separate out gradient volumes
      this->m_IntensityData = DWIVectorImageType::New();
      this->m_IntensityData->SetNumberOfComponentsPerPixel(gradientIndices.size());
      this->m_IntensityData->CopyInformation(originalNrrd);
      this->m_IntensityData->SetRegions(region);
      this->m_IntensityData->Allocate();

      const unsigned int numGradientIndices = gradientIndices.size();
      //
      // make new vectorimage with just the gradient values
      itk::ImageRegionIterator<DWIVectorImageType> newIt(this->m_IntensityData,region);
      for(b0It.Begin(), newIt.Begin(); !b0It.IsAtEnd() && !newIt.IsAtEnd(); ++b0It,++newIt)
        {
        DWIVectorImageType::PixelType newVal;
        const DWIVectorImageType::PixelType oldVal = b0It.Get();
        for(unsigned int i = 0; i < numGradientIndices; ++i)
          {
          newVal[i] = oldVal[gradientIndices[i]];
          }
        newIt.Set(newVal);
        }

      // get just the non-B0 gradients
      for(unsigned int i = 0; i < numGradientIndices; ++i)
        {
        this->m_VoxelLatticeAlignedGradientDirections.conservativeResize(i+1,3);
        this->m_VoxelLatticeAlignedGradientDirections.row(i) = this->m_Gradients.row(gradientIndices[i]);
        }
      return true;
    }
  bool Compute()
    {
      if(!this->AverageB0AndExtractIntensity())
        {
        return false;
        }
      // MeasurementFrame was a parameter in the original computation
      // but isn't needed here -- gradients are converted to world
      // frame when the file gets loaded.

      unsigned numGradientDirections = this->m_IntensityData->GetNumberOfComponentsPerPixel();

      //
      // MATLAB CODE IS
      // DWIIntensityData = single(DWIIntensityData) ./ averagedB0(:,:,:,ones(1,numGradientDirections));
      // %remove negative values
      // DWIIntensityData(DWIIntensityData<0)=eps;
      itk::ImageRegionIterator<DWIVectorImageType>
        ivIt(this->m_IntensityData,this->m_IntensityData->GetLargestPossibleRegion());
      unsigned int numGrad = this->m_IntensityData->GetNumberOfComponentsPerPixel();

      itk::ImageRegionConstIterator<B0AvgImageType>
        b0It(this->m_AverageB0,this->m_AverageB0->GetLargestPossibleRegion());
      for(ivIt.Begin(),b0It.Begin(); !ivIt.IsAtEnd() && !b0It.IsAtEnd(); ++ivIt,++b0It)
        {
        DWIVectorImageType::PixelType vec = ivIt.Get();
        double curB0 = b0It.Get();
          if(std::fabs(curB0) > eps)
            {
            for(unsigned int i = 0; i < numGrad; ++i)
              {
              vec[i] /= curB0;
              if(vec[i] < 0)
                {
                vec[i] = eps;
                }
              }
            }
          else
            {
            for(unsigned int i = 0; i < numGrad; ++i)
              {
              vec[i] = eps;
              }
            }
          ivIt.Set(vec);
        }

      MatrixType estimatedGradients = this->m_Gradients;
      // this matlab code doesn't seem to do anything:
      // n0 = size(new_gradients,1);
      // new_gradients=new_gradients(1:n0,:);
      DWIVectorImageType::SizeType dWIVecImageSize =
        this->M_IntensityData->GetLargestPossibleRegion().GetSize();
      unsigned numGradientVoxels =
        dwiVecImageSize[0]
        * dwiVecImageSize[1]
        * dwiVecImageSize[2];
      MatrixType icos3;
      Icosahedron3(icos3);
      MatrixType D0(3,3) = MatrixType::Zeros(3,3);
      
      double rho;
      double phi;
      this->OptimBand(m_BValue,D0,icos3,rho,phi);
      return true;
    }
  void OptimBand(double B,const MatrixType &D,const MatrixType &ico3,double &rho, double &phi)
    {
      MatrixType u = ico3;

      // S=exp(-b*sum((u*D).*u,2));
      MatrixType u2 = u*D;
      MatrixType u3(u.rows(),u.cols());

      MatrixType S(u.rows(),1);

      for(unsigned i = 0; i < u.rows(); ++i)
        {
        for(unsigned j = 0; j < u.cols(); ++j)
          {
          u3(i,j) = u2(i,j) * u(i,j);
          }
        }
      for(unsigned i = 0; i < u.rows(); ++i)
        {
        S(i) = std::exp(-b * u3.row(i).sum());
        }
      MatrixType Y = SphericalHarmonics(u,22);
      
    }
private:
  NrrdFile &                    m_NrrdFile;
  MaskImageType::Pointer        m_MaskImage;
  DWIVectorImageType::PointType m_SpaceOrigin;
  MatrixType                    m_Gradients;
  double                        m_BValue;
  MatrixType                    m_SpaceDirections;
  MatrixType                    m_VoxelLatticeAlignedGradientDirections;
  DWIVectorImageType::Pointer   m_IntensityData;
  B0AvgImageType::Pointer       m_AverageB0;
};
#endif // __DoCSEstimate_h
