#ifndef __BalancedDWIReplications_h
#define __BalancedDWIReplications_h
#include "nrrdIO.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include <cmath>
#include <algorithm>
class BalancedDWIReplications
{
public:
  BalancedDWIReplications(NrrdFile &rawDWI) : m_NrrdFile(rawDWI) {}

  void ReplicateGradient(unsigned gdindex, unsigned count)
    {

      MatrixType newGradients(this->m_NrrdFile.GetGradients());

      // duplicate gradients
      for(unsigned int i = 0; i < count; ++i)
        {
        unsigned rows = newGradients.rows();
        newGradients.conservativeResize(rows+1,3);
        newGradients.row(rows) = newGradients.row(gdindex);
        }

      this->m_NrrdFile.SetGradients(newGradients);;

      //
      // duplicate image data, using itk::PasteImageFilter
      NrrdFile::ImageType::Pointer oldImage = m_NrrdFile.GetImage();
      unsigned int oldVectorSize = oldImage->GetNumberOfComponentsPerPixel();
      unsigned int newVectorSize = oldVectorSize+count;

      NrrdFile::ImageType::Pointer newImage = NrrdFile::ImageType::New();
      newImage->SetSpacing(oldImage->GetSpacing());
      newImage->SetDirection(oldImage->GetDirection());
      newImage->SetOrigin(oldImage->GetOrigin());
      newImage->SetRegions(oldImage->GetLargestPossibleRegion());
      newImage->SetNumberOfComponentsPerPixel(newVectorSize);
      newImage->Allocate();

      typedef itk::ImageRegionIterator<NrrdFile::ImageType> ImageIterator;
      typedef itk::ImageRegionConstIterator<NrrdFile::ImageType> ImageConstIterator;

      ImageConstIterator sourceIt(oldImage,oldImage->GetLargestPossibleRegion());
      ImageIterator      destIt(newImage,newImage->GetLargestPossibleRegion());

      for(sourceIt.GoToBegin(), destIt.GoToBegin(); !sourceIt.IsAtEnd(); ++sourceIt, ++destIt)
        {
        const NrrdFile::ImageType::PixelType oldVal = sourceIt.Get();
        NrrdFile::ImageType::PixelType newVal;
        unsigned int i = 0;
        for(;i < oldVectorSize; ++i)
          {
          newVal[i] = oldVal[i];
          }
        for(;i < newVectorSize; ++i)
          {
          newVal[i] = oldVal[gdindex];
          }
        destIt.Set(newVal);
        }

      unsigned rows = this->m_Metric.rows();
      this->m_Metric.conservativeResize(rows+count);
        // update the metrics array
      for(unsigned int i = rows; i < (rows+count); ++i)
        {
        this->m_Metric.row(i) = this->m_Metric.row(gdindex);
        }

      this->m_NrrdFile.SetImage(newImage);
    }
  void compute()
    {
      typedef Eigen::Array<unsigned,Eigen::Dynamic,1> ArrayType;

      const double Pi(std::atan(static_cast<double>(1.0))*4);

      const MatrixType gradients = m_NrrdFile.GetGradients();
      unsigned int n = gradients.rows();

      this->m_Metric.conservativeResize(n);
      this->m_Metric.setConstant(90);

      ArrayType isGradient(n,1);
      isGradient.setOnes();

      ArrayType isLess5Degrees(n,1);
      isLess5Degrees.setZero();

      ArrayType counts(n,1);
      counts.setOnes();

      for(unsigned int i = 0; i < n; ++i)
        {
        double norm = gradients.row(i).norm();
        if(norm < 0.00001)
          {
          isGradient(i) = 0;
          isLess5Degrees(i) = 0;
          continue;
          }
        for(unsigned j = 0; j < n; ++j)
          {
          if(i != j)
            {
            Eigen::Vector3d rowI = gradients.row(i);
            Eigen::Vector3d rowJ = gradients.row(j);
            double dot = std::abs(rowI.dot(rowJ));
            double angleBetween = (std::acos(dot) * 180) / Pi;
            this->m_Metric(i) = std::min(this->m_Metric(i),angleBetween);
            this->m_Metric(j) = std::min(this->m_Metric(j),angleBetween);
            if(angleBetween < 5)
              {
              isLess5Degrees(i) = 1;
              isLess5Degrees(j) = 1;
              counts(i) = counts(i) + 1;
              counts(j) = counts(j) + 1;
              }
            }
          }
        }
      unsigned _max = counts.maxCoeff();
      // TODO this expression is insane. It isn't really possible to
      // decode, though I try.
      // GradientsNeedDuplications =
      // ( (counts < max(counts) ).*max(counts) - counts .* (counts < max(counts) ) ) .* isGradient;
      ArrayType GradientsNeedDuplications(n,1);
      for(unsigned int i = 0; i < counts.cols(); ++i)
        {
        unsigned tmp = 0;
        if(counts(i) < _max)
          {
          tmp = _max;
          }
        unsigned int tmp2 = 0;
        if(counts(i) > 0 && counts(i) <  _max)
          {
          tmp2 = counts(i);
          }
        GradientsNeedDuplications(i) = (tmp - tmp2) * isGradient(i);
        }
      for(unsigned i = 0; i < n; ++i)
        {
        if(GradientsNeedDuplications(i) > 0)
          {
          this->ReplicateGradient(i,GradientsNeedDuplications(i));
          }
        }
    }
  const VectorType &GetMetric() const { return m_Metric; }
private:
  NrrdFile &m_NrrdFile;
  VectorType m_Metric;
};

#endif // __BalancedDWIReplications_h
