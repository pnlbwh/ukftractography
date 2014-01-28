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

class NrrdFile
{
public:
  typedef itk::VectorImage<double,3>      ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef std::vector<double>             BValueVector;
  NrrdFile() : m_GradientCount(0), m_BValue(-1.0)
    {}

  void
  ExtractGradients()
    {
      const itk::MetaDataDictionary &thisDic = this->m_Image->GetMetaDataDictionary();
      // determine if image contains bvalue & gradients.
      std::string bValueString;

      if(itk::ExposeMetaData<std::string>(thisDic,"DWMRI_b-value",bValueString) == true)
        {
        std::stringstream ss(bValueString);
        ss >> this->m_BValue;

        std::stringstream gradTag;
        m_GradientCount = 0;
        gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << this->m_GradientCount;

        std::string gradString;

        while(itk::ExposeMetaData<std::string>(thisDic,gradTag.str(),gradString) != false)
          {
          VectorType curGrad(3,1);
          double val;

          std::stringstream grad(gradString);

          for(unsigned i = 0; i < 3; ++i)
            {
            grad >> val;
            curGrad(i) = val;
            }
          this->m_Gradients.col(this->m_GradientCount) = curGrad;
          gradTag.str("");
          ++this->m_GradientCount;
          gradTag << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << this->m_GradientCount;
          }
        }
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
      ExtractGradients();
      return true;
    }

  void
  SaveFile(const std::string &filename)
    {
      itk::MetaDataDictionary & thisDic = this->m_Image->GetMetaDataDictionary();
      //
      // space direction
      itk::EncapsulateMetaData<std::string>(thisDic,"NRRD_space","left-posterior-superior");
      itk::EncapsulateMetaData<std::string>(thisDic,std::string("modality"),
                                            std::string("DWMRI"));

      std::stringstream ss;
      ss.precision(10);
      ss << this->m_BValue;
      itk::EncapsulateMetaData<std::string>(thisDic,std::string("DWMRI_b-value"),
                                            ss.str());

      for(unsigned i = 0; i < this->m_Gradients.cols(); ++i)
        {
        ss.str(std::string());
        ss << std::setprecision(17) << std::scientific
           << this->m_Gradients(0,i) << " "
           << this->m_Gradients(1,i) << " "
           << this->m_Gradients(2,i);
        std::stringstream gradName;
        gradName << "DWIMRI_gradient_" << std::setw(4) << std::setfill('0') << i;
        itk::EncapsulateMetaData<std::string>(thisDic,
                                              gradName.str(),ss.str());
        }

      WriterType::Pointer writer = WriterType::New();
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
private:
  ImageType::Pointer m_Image;
  MatrixType m_Gradients;
  unsigned int m_GradientCount;
  double m_BValue;
};
#endif // __nrrdIO_h
