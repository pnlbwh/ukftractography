#ifndef VTKREADER_H_
#define VTKREADER_H_

#include <string>
#include <vector>
#include <map>
#include "linalg.h"
#include "CompareFiber.h"
#include <fstream>

// Class that allows to read a .vtk file.
class vtkReader
{
public:
  vtkReader() : _sInputPath(""), _fibers(NULL)
  {
  }

  virtual ~vtkReader()
  {
  }

  bool Run();

  void SetInputPath(const std::string & path);

  void SetOutputFibers(std::vector<Fiber> & fibers);

  void SetReadFieldData(const bool option);

protected:
private:
  std::string          _sInputPath;
  std::vector<Fiber> * _fibers;

  // internal parameters
  std::vector<std::vector<unsigned int> > _lines;

  int  _nNumOfPoints;
  int  _nNumOfFibers;
  int  _nNumOfFields;
  bool _bReadFieldData;

  bool ReadLines(std::ifstream & input);

  bool ReadRest(std::ifstream & input);

};

#endif  // VTKREADER_H_
