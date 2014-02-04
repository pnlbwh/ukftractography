#ifndef __nrrdIO_h
#define __nrrdIO_h
#include <fstream>
#include "CompressedSensing.h"
#include "itkMetaDataObject.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <ctype.h>

/** NrrdFile holds original volume + gradients
 *
 */
class NrrdFile
{
public:
  typedef signed short                    PixelType;
  typedef itk::VectorImage<PixelType,3>   ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef std::vector<double>             BValueVector;

  NrrdFile() : m_GradientCount(0), m_BValue(-1.0)
    {}

  bool
  ExtractGradients()
    {

      itk::MetaDataDictionary &thisDic = this->m_Image->GetMetaDataDictionary();
      // determine if image contains bvalue & gradients.
      std::string bValueString;

      if(itk::ExposeMetaData<std::string>(thisDic,"DWMRI_b-value",bValueString) == false)
        {
        std::cerr << "Gradient information missing from this image" << std::endl;
        return false;
        }

      std::stringstream ss(bValueString);
      ss >> this->m_BValue;

      std::stringstream gradTag;
      m_GradientCount = 0;
      gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << this->m_GradientCount;

      std::string gradString;

      while(itk::ExposeMetaData<std::string>(thisDic,gradTag.str(),gradString) != false)
        {
        this->m_Gradients.conservativeResize(this->m_GradientCount+1,3);
        VectorType curGrad(3,1);
        double val;

        std::stringstream grad(gradString);

        for(unsigned i = 0; i < 3; ++i)
          {
          grad >> val;
          curGrad(i) = val;
          }
        this->m_Gradients.row(this->m_GradientCount) = curGrad;
        gradTag.str("");
        ++this->m_GradientCount;
        gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << this->m_GradientCount;
        }
      //
      // remove the measurment frame
      Eigen::Matrix3d mf;
      std::vector<std::vector<double> > msrFrame;
      if(itk::ExposeMetaData<std::vector<std::vector<double> > >( thisDic, "NRRD_measurement frame", msrFrame ) == false)
        {
        return true;
        }
      for(unsigned int i = 0; i < 3; ++i)
        {
        for(unsigned int j = 0; j < 3; ++j)
          {
          mf(i,j) = msrFrame[i][j];
          }
        }
      mf = mf.inverse().eval();
      for(unsigned i = 0; i < this->m_GradientCount; ++i)
        {
        Eigen::Vector3d curRow = this->m_Gradients.row(i);
        curRow = mf * curRow;
        if(curRow.norm() > 0.0)
          {
          curRow.normalize();
          }
        this->m_Gradients.row(i) = curRow;
        }
      for(unsigned int i = 0; i < 3; ++i)
        {
        for(unsigned int j = 0; j < 3; ++j)
          {
          msrFrame[i][j] = i == j ? 1.0 : 0.0;
          }
        }

      itk::EncapsulateMetaData<std::vector<std::vector<double> > >(thisDic, "NRRD_measurement frame", msrFrame );
      return true;
    }

  bool LoadFile(const std::string &filename)
    {
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(filename);
      try
        {
        reader->Update();
        }
      catch(itk::ExceptionObject &excp)
        {
        std::string msg("itkLoadWithMetaData: can't read ");
        msg += filename;
        msg += " ";
        msg += excp.what();
        std::cerr << msg;
        return false;
        }
      this->m_Image = reader->GetOutput();
      if(!ExtractGradients())
        {
        return false;
        }
      return true;
    }

  void
  SaveFile(const std::string &filename)
    {
      itk::MetaDataDictionary & thisDic = this->m_Image->GetMetaDataDictionary();
      //
      // space direction
      itk::EncapsulateMetaData<std::string>(thisDic,"NRRD_space","left-posterior-superior");
      itk::EncapsulateMetaData<std::string>(thisDic,std::string("modality"),std::string("DWMRI"));

      std::stringstream ss;
      ss.precision(10);
      ss << this->m_BValue;
      itk::EncapsulateMetaData<std::string>(thisDic,std::string("DWMRI_b-value"),
                                            ss.str());

      for(unsigned i = 0; i < this->m_Gradients.rows(); ++i)
        {
        ss.str(std::string());
        ss << std::setprecision(17) << std::scientific
           << this->m_Gradients(i,0) << " "
           << this->m_Gradients(i,1) << " "
           << this->m_Gradients(i,2);
        std::stringstream gradName;
        gradName << "DWIMRI_gradient_" << std::setw(4) << std::setfill('0') << i;
        itk::EncapsulateMetaData<std::string>(thisDic,gradName.str(),ss.str());
        }

      WriterType::Pointer writer = WriterType::New();
      writer->UseCompressionOn();
      writer->UseInputMetaDataDictionaryOn();
      writer->SetFileName(filename);
      writer->SetInput(this->m_Image);
      try
        {
        writer->Write();
        }
      catch(itk::ExceptionObject &excp)
        {
        std::string msg("itkSaveWithMetaData: can't read ");
        msg += filename;
        msg += " ";
        msg += excp.what();
        std::cerr << msg << std::endl;
        }
    }
  MatrixType GetGradients() const { return m_Gradients; }
  void SetGradients(const MatrixType &newGrad) { m_Gradients = newGrad; }
  ImageType::Pointer GetImage() { return m_Image; }
  void SetImage(ImageType::Pointer &newImage) { m_Image = newImage; }
  double GetBValue() { return m_BValue; }
  void SetBValue(double newB) { m_BValue = newB; }
private:
  ImageType::Pointer m_Image;
  MatrixType m_Gradients;
  unsigned int m_GradientCount;
  double m_BValue;
};
#endif // __nrrdIO_h
