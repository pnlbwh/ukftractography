/**
 * \file Region.h
 * \brief Contains class Region which defines a ROI that can be used for logical fiber expressions.
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef REGION_H_
#define REGION_H_

#include <string>
#include <vector>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>

/**
 * \class Region
 * \brief Contains definition of a ROI of interest for filtering out fibers
 */
class Region
{
public:

  /** Defines type of a 3D unsigned integer image */
  typedef itk::Image< unsigned int, 3> UIntImageType;

  /** 
   * Constructor 
   * \param path Path of a label file
   * \param labels Labels of interest of that file
  */
  Region(const std::string & path, const std::vector<int> & labels);
  
  /** Virtual Destructor */
  virtual ~Region() { };

  /** Return X Dimension of the Region */
  int SizeX() const {
    return _nSizeX;
  }

  /** Return Y Dimension of the Region */
  int SizeY() const {
    return _nSizeY;
  }

  /** Return Z Dimension of the Region */
  int SizeZ() const {
    return _nSizeZ;
  }

  /** Return the Look-up table */
  const std::vector<int>& GetLut() const {
    return _lut;
  }

  /** Return if the region is empty, i.e. no path to label file is given */
  const bool IsEmpty() const {
    return _bEmpty;
  }

  /** Reference image, used for reading with ITK */
  static UIntImageType::Pointer refImage;

protected:

  /** checks wheather the integer n is in the vector v */
  bool inVector(const int & n, const std::vector<int> & v);

  /** X Dimension of Region */
  int _nSizeX;

  /** Y Dimension of Region */
  int _nSizeY; 

  /** Z Dimension of Region */
  int _nSizeZ;

  /** A look-up table for quick access to region elements */
  std::vector<int> _lut;

  /** Wheater the label path is empty */
  bool _bEmpty;

};

#endif  // REGION_H_
