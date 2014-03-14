#ifndef __NrrdFile_h
#define __NrrdFile_h
#include "CompressedSensing.h"
#include "itkVectorImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkMetaDataObject.h"
#include <vector>
/** NrrdFile holds original volume + gradients
 *
 */
class NrrdFile
{
public:
  typedef double                          PixelType;
  typedef itk::VectorImage<PixelType,3>   ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typedef std::vector<double>             BValueVector;

  NrrdFile() : m_GradientCount(0), m_BValue(-1.0)
    {}

  bool ExtractGradients();
  bool LoadFile(const std::string &filename);
  void SaveFile(const std::string &filename);
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
  itk::MetaDataDictionary m_SaveDict;
};
#endif // __NrrdFile_h
