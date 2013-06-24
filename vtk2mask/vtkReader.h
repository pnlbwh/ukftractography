/**
 * \file vtkReader.h
 * \brief Contains class definition of vtkReader
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com) based on code by Stefan Lienhard
*/

#ifndef VTKREADER_H_
#define VTKREADER_H_

#include <string>
#include <vector>
#include <map>
#include "linalg.h"
#include "fiber.h"
#include <fstream>

/**
 * \class vtkReader
 * \brief Class that allows to read a .vtk file.
*/
class vtkReader
{
public:
  /** Constructor */
  vtkReader()
  {
  }

  /** Destructor */
  virtual ~vtkReader()
  {
  }

  /** Start the reading of the file with the specified values */
  bool Run();

  /** Set the path of the vtk file to be read */
  void SetInputPath(const std::string & path);

  /** Set the pointer to the fiber vector to store the data from the vtk file */
  void SetOutputFibers(std::vector<Fiber> & fibers);

  /** Set whether the additional fields should be read or not */
  void SetReadFieldData(const bool option);

  /** Set wheater the program should outputput execution infos on the commandline */
  void SetVerbose(bool b)
  {
    _bVerbose = b;
  }

protected:

  /** Pointer to the input path of the vtk File */
  const char * _sInputPath;

  /** Pointer to the fiber vector where the input is stored */
  std::vector<Fiber> * _fibers;

  /** Array to temporalily store the lines field of the vtk */
  std::vector<std::vector<unsigned int> > _lines;

  /** Total number of points in the vtk */
  int _nNumOfPoints;

  /** Total number of fibers in the vtk */
  int _nNumOfFibers;

  /** Total number of scalar fields in the vtk */
  int _nNumOfFields;

  /** Wheater to read the field data */
  bool _bReadFieldData;

  /** Wheater to output execution infos on the command line */
  bool _bVerbose;

  /** Read the Line information of the vtk, needs to be done before reading points. */
  bool ReadLines(std::ifstream & input);

  /** Read the points data, and field data */
  bool ReadRest(std::ifstream & input);

};

#endif  // VTKREADER_H_
