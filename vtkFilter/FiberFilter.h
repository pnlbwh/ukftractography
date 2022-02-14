/**
 * \file FiberFilter.h
 * \brief Contains class FiberFilter for filtering VTK fibers that go through a certain Region of interest
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef FIBERFILTER_H_
#define FIBERFILTER_H_

#include <string>
#include <vector>
#include <map>
#include "linalg.h"
#include "fiber.h"
#include "Region.h"

#include <itkImage.h>
#include <itkMacro.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>

/**
 * \class FiberFilter
 * \brief Provides the functinality for actually filtering fibers
*/
class FiberFilter
{
public:

  /** Unisgned image 3D Image type */
  typedef itk::Image<unsigned int, 3> UIntImageType;

  /** There are modes, 'connecting to', or 'ending in', where there is no distinction between 'ending in', and 'starting
    in' */
  enum ConnectionMode
    {
    PASS = true,
    END = false
    };

  /** Constructor */
  FiberFilter();

  /** Virtual Destructor */
  virtual ~FiberFilter()
  {
  }

  /** Run the fiber filter when all parameters are set */
  bool Run();

  /** Set the pointer to the input fiber vector */
  void SetInputFibers(const std::vector<Fiber> & fibers);

  /** Set the pointer to the output fiber vector */
  void SetOutputFibers(std::vector<Fiber> & fibers);

  /** Set the pointer to a region of interest */
  void SetRegion(const Region & region);

  /** Set whether to run in 'connecting to', or 'ending in' mode */
  void SetConnectionMode(bool option);

  /** Whether to filter the given fibers in inverse mode. I.e. output all fibers that are n o t in the region */
  void SetCalcInverse(bool option);

protected:

  /** Pointer to the region look-up table */
  const int* _lut;

  /** Pointer to a reference image */
  UIntImageType::Pointer _refImage;

  /** Size of region in X */
  int _nx;

  /** Size of region in X */
  int _ny;

  /** Size of region in X */
  int _nz;

  /** Size of the look-up table */
  unsigned int _nSizeLUT;

  /** Pointer to the input fibers */
  const std::vector<Fiber> * _inFibers;

  /** Pointer to the output fibers */
  std::vector<Fiber> * _outFibers;

  /** Poiter to the region which should be filtered */
  const Region * _region;

  /** Were the input fibers set? */
  bool _bInSet;

  /** Were the output fibers set? */
  bool _bOutSet;

  /** Was the region set? */
  bool _bRegionSet;

  /** connecting or ending in mode */
  bool _bConnectionMode;

  /** Inverse, or normal mode */
  bool _bCalcInverse;

  /** Check if the fiber is in the region */
  bool inRegion(const Fiber & fiber);

  /** returns inRegion in normal mode, and !inRegion for inverse mode */
  bool CheckConditions(const Fiber & fiber);

};

#endif  // FIBERFILTER_H_
