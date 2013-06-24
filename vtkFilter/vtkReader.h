#ifndef VTKREADER_H_
#define VTKREADER_H_

#include <string>
#include <vector>
#include <map>
#include "linalg.h"
#include "fiber.h"
#include <fstream>

class vtkReader
{
public:
  vtkReader()
  {
  }

  virtual ~vtkReader()
  {
  }

  bool Run();

  void SetInputPath(const std::string & path);

  void SetOutputFibers(std::vector<Fiber> & fibers);

  void SetReadFieldData(const bool option);

  void SetVerbose(bool b)
  {
    _bVerbose = b;
  }

protected:

  const char *         _sInputPath;
  std::vector<Fiber> * _fibers;

  // internal parameters
  std::vector<std::vector<unsigned int> > _lines;

  int  _nNumOfPoints;
  int  _nNumOfFibers;
  int  _nNumOfFields;
  bool _bReadFieldData;
  bool _bVerbose;

  bool ReadLines(std::ifstream & input);

  bool ReadRest(std::ifstream & input);

};

#endif  // VTKREADER_H_
