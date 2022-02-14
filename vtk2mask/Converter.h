/**
 * \file Converter.h
 * \brief Contains functionality to go from a VTK fiber to a NRRD Mask
*/

#ifndef CONVERTER_H_
#define CONVERTER_H_

#include <cassert>
#include <string>
#include "Array3D.h"
#include "fiber.h"

#include <itkImage.h>
#include <itkMacro.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>

/**
 * \class Converter
 * \brief Functionality for converting a fiber vector to a scalar mask
 *
 * Converts each point of the fibers to the reference dwi space and
 * stores a scalar value of the fiber at the respective map point
*/
class Converter
{
public:

  /** A 3D ITK float vector, basic component of the dwi image */
  typedef itk::Vector<float, 3> VectorType;

  /** Type of a float diffusion weighted MR Image */
  typedef itk::Image<VectorType, 3> DiffusionImageType;

  /** Type of a scalar float image */
  typedef itk::Image<float, 3> FloatImageType;

  /** Type of a scalar unsigned int image */
  typedef itk::Image<unsigned int, 3> UIntImageType;

  /** A single float point */
  typedef itk::Point<float, FloatImageType::ImageDimension> FloatPointType;

  /** Constructor */
  Converter();

  /** Virtual destructor */
  virtual ~Converter()
  {
  }

  /** Run the conversion when all parameters are set */
  bool Run();

  /** Set the fibers to be converted */
  void SetInputFibers(std::vector<Fiber> & fibers);

  /** Set the path to the reference DWI Volume, which is used to define the space of the new scalar map */
  void SetReferenceFile(std::string & reference_volume);

  /** Set the path to the Output NRRD volume */
  void SetOutputVolumeFile(std::string & output_volume);

  /** Set the path to the Standart Deviation Output NRRD volume */
  void SetStdDevFile(std::string & std_dev_volume);

  /** Set wheather to ouput execution informations on the command line */
  void SetVerbose(bool verbose);

  /** Set the label file for only performing the conversion in certain areas */
  void SetLabelFile(std::string & label_file);

  /** Select Labels from the label file */
  void SetLabelNumber(int label);

  /** Set what field should be extracted from the fiber, or binary map if left empty */
  void SetFieldName(std::string & field_name);

protected:

  /** Pointer to the input fiber vector */
  std::vector<Fiber> * _fibers;

  /** Output file path of the scalar map */
  FloatImageType::Pointer _nrrdDataOut;

  /** Output file path of the standart deviation map */
  FloatImageType::Pointer _nrrdVarDataOut; // only allocate if Variance are required

  /** A temp image for the label information */
  UIntImageType::Pointer _Label;

  /** A large 3D Array where each point is a vector (so really 4D). Each time a fiber passes through the point its
    scalar field value is pushed_back to the vector at this point. */
  Array3D<std::vector<float> > * _matField;

  /** The dimensions of the data in X direction */
  int _nDimX;

  /** The dimensions of the data in Y direction */
  int _nDimY;

  /** The dimensions of the data in Z direction */
  int _nDimZ;

  /** Wheather to output the Standart Deviation map */
  bool _bOutputVarianceMaps;

  /** Wheater a scalar field was selected from the VTK, if false a scalar map of where the fibers pass through will be
    generated. */
  bool _bOutputField;

  /** Wheater a label file was given */
  bool _bHasLabel;

  /** Output execution infos on the command line */
  bool _bVerbose;

  /** Path of the Output Scalar NRRD Volume */
  std::string _sOutputVolume;

  /** Path of the Reference DWI Volume */
  std::string _sReferenceVolume;

  /** Path of the Output Standart Deviation Volume */
  std::string _sStandartDevVolume;

  /** Path of the Label File */
  std::string _sLabelFile;

  /** Name of the Scalar Field of interest */
  std::string _sScalarFieldName;

  /** Which label of the label file are we interested in */
  int _nLabelOfInterest;

  /** Wheater to output the mean, and standart deviation of the field */
  void OutputStats();

  /** Allocate the output NRRDs */
  void CreateOutNrrd(const std::string & nrrd_path);

  /** Write data to the Output NRRD */
  void WriteOutNrrd(const std::string & out_path);

  /**
   * Fill a large 3D Matrix with the fiber scalar values at each where the
   * indices of the matrix correspond to the ijk space.
   * \param field_name Which field to process
  */
  bool FillMatField(const std::string & field_name);

  /** Basically flaten the _matField, i.e. average values at each point */
  void AverageVoxels();

  /** Load the Label File */
  void LoadLabel();

  /**
  * Extract the scalar values from the Fiber, average each point, and write the output
  * \param field_name Which field to process
  * \param path Path of the output NRRD File.
  */
  void ProcessField(const std::string & field_name, const std::string & path);

  /**
   * Calculate the standart deviation of a vector
   * \param vec The vector
   * \param mean The mean of the vector
  */
  float CalcStdDev(const std::vector<float> & vec, const float & mean);

  /** Calculate the mean of a vector */
  float CalcMean(const std::vector<float> & vec);

};

#endif  // CONVERTER_H_
