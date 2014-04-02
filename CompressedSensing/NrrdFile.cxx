#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <ctype.h>
#include "NrrdFile.h"

bool
NrrdFile
::ExtractGradients()
{

  this->m_SaveDict  = this->m_Image->GetMetaDataDictionary();
  // determine if image contains bvalue & gradients.
  std::string bValueString;

  if(itk::ExposeMetaData<std::string>(this->m_SaveDict,"DWMRI_b-value",bValueString) == false)
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

  while(itk::ExposeMetaData<std::string>(this->m_SaveDict,gradTag.str(),gradString) != false)
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
  Eigen::Matrix3d mf;
  std::vector<std::vector<double> > msrFrame;
  if(itk::ExposeMetaData<std::vector<std::vector<double> > >( this->m_SaveDict, "NRRD_measurement frame", msrFrame ) == false)
    {
    itk::EncapsulateMetaData<std::vector<std::vector<double> > >(this->m_SaveDict, "NRRD_measurement frame", msrFrame );
    return true;
    }
  ImageType::DirectionType directions = this->m_Image->GetDirection();
  Eigen::Matrix3d dirCosines;

  for(unsigned int i = 0; i < 3; ++i)
    {
    for(unsigned int j = 0; j < 3; ++j)
      {
      mf(i,j) = msrFrame[i][j];
      dirCosines(i,j) = directions[i][j];
      }
    }
  // HANS TODO -- what frame of reference should the gradients
  // live in? The #if 0 code is what I cribbed out of DWIConvert,
  // which  removes the measurement frame from the gradients.
  // what the CompressedDWI Matlab code does is multiply each
  // gradient by the inverse of the dir cosines and the
  // measurement frame.  That's what I do below in order to be
  // consistent with the matlab code
#if 0
  //
  // remove the measurment frame
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
#else
  dirCosines = dirCosines.inverse().eval();
  mf = dirCosines * mf;
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
#endif
  itk::EncapsulateMetaData<std::vector<std::vector<double> > >(this->m_SaveDict, "NRRD_measurement frame", msrFrame );
  return true;
}
bool
NrrdFile
::LoadFile(const std::string &filename)
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
NrrdFile
::SaveFile(const std::string &filename)
{

  itk::MetaDataDictionary & thisDic = this->m_Image->GetMetaDataDictionary();
  //
  // copy dictionary except for the gradients
  for(itk::MetaDataDictionary::ConstIterator it = this->m_SaveDict.Begin();
      it != this->m_SaveDict.End(); ++it)
    {
    if(it->first.find("DWMRI_gradient_",0) == std::string::npos)
      {
      thisDic.Set(it->first,it->second);
      }
    }

  //
  // space direction
  itk::EncapsulateMetaData<std::string>(thisDic,"NRRD_space","left-posterior-superior");

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
