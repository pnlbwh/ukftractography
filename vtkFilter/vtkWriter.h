/**
 * \file vtkWriter.h
 * \brief Contains class vtkWriter for Writing a fiber vector to a vtk file
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef VTKWRITER_H_
#define VTKWRITER_H_

#include <string>
#include <vector>
#include <map>
#include "linalg.h"
#include "fiber.h"

/**
 * \class vtkWriter
 * \brief Contains functions to write a vector of Fibers to a *.vtk file
*/
class vtkWriter
{
public:

  /** Constructor */
  vtkWriter()
  {
  }

  /** Virtual Destructor */
  virtual ~vtkWriter()
  {
  }

  /** Run the writer, when all parameters are set */
  bool Run();

  /** Set the Output vtk path */
  void SetOutputPath(const std::string& path);

  /** Set the pointer to the fiber vector */
  void SetInputFibers(std::vector<Fiber> & fibers);

protected:

  /** The output path of the vtk file */
  const char * _sOutputPath;

  /** Pointer to the fiber vector */
  std::vector<Fiber> * _fibers;

  /** Number of points of all fibers */
  size_t _nNumOfPoints;

  /** Number of fibers */
  size_t _nNumOfFibers;

  /** Number of scalar fields of the fiber */
  size_t _nNumOfFields;

  /** The respective lengths of each fiber */
  std::vector<size_t> _fiberLengths;

  /** Write the top line of the vtk file */
  void WriteHeader(std::ofstream & output);

  /** Write the points section */
  void WritePoints(std::ofstream & output);

  /** Write the lines section */
  void WriteLines(std::ofstream & output);

  /** Write the scalar fields */
  void WriteFields(std::ofstream & output);

};

#endif  // VTKWRITER_H_
