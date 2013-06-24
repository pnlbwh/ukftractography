#include "vtkReader.h"
#include <fstream>
#include <iostream>

bool vtkReader::Run()
{

  std::ifstream input(_sInputPath.c_str() );

  if( !input.is_open() )
    {
    std::cout << "Failed to open " << _sInputPath << "." << std::endl;
    return 1;
    }

  if( ReadLines(input) )
    {
    std::cout << "Couldn't read LINES\n";
    return 1;
    }

  if( ReadRest(input) )
    {
    std::cout << "Error parsing the vtk file.\n";
    return 1;
    }
  input.close();
  return 0;
}

void vtkReader::SetInputPath(const std::string & path)
{
  _sInputPath = path;
}

void vtkReader::SetOutputFibers(std::vector<Fiber> & fibers)
{
  _fibers = &fibers;
}

void vtkReader::SetReadFieldData(const bool option)
{
  _bReadFieldData = option;
}

bool vtkReader::ReadLines(std::ifstream & input)
{

  input.seekg(std::ios::beg);

  std::string sCurr;
  int         nCurr;

  while( !input.eof() )
    {

    // skip forward to lines
    while( sCurr != "LINES" )
      {
      input >> sCurr;
      }

    input >> nCurr;
    _nNumOfFibers = nCurr;
    _lines.resize(_nNumOfFibers);
    _fibers->resize(_nNumOfFibers);

    input >> nCurr;
    _nNumOfPoints = nCurr - _nNumOfFibers;

    int nLineLength;
    for( int i = 0; i < _nNumOfFibers; ++i )
      {
      // first element of line is fiber length
      input >> nLineLength;
      _lines[i].resize(nLineLength);
      // the remaining elements are the point indices
      for( int j = 0; j < nLineLength; ++j )
        {
        input >> _lines[i][j];
        }
      }
    return 0;
    }

  return 1;
}

bool vtkReader::ReadRest(std::ifstream & input)
{

  // TODO: actually use _lines data because vtk fibers might be unordered
  input.seekg(std::ios::beg);

  std::string sCurr;
  int         nCurr;
  int         num_points = 0;

  while( !input.eof() )
    {
    while( input.peek() == '#' )
      {
      getline(input, sCurr);
      }

    input >> sCurr;
    // std::cout << sCurr << " ";

    if( !sCurr.compare("DATASET") )
      {
      input >> sCurr;
      assert(!sCurr.compare("POLYDATA") );
      input >> sCurr;
      assert(!sCurr.compare("POINTS") );
      std::cout << "Found Point data...\n";

      input >> num_points;
      assert(num_points == _nNumOfPoints);

      input >> sCurr;
      assert(!sCurr.compare("float") );

      float x, y, z;
      for( int i = 0; i < _nNumOfFibers; ++i )
        {
        int nLineLength = _lines[i].size();
        (*_fibers)[i].Points.resize(nLineLength);
        for( int j = 0; j < nLineLength; ++j )
          {
          input >> x;
          input >> y;
          input >> z;
          (*_fibers)[i].Points[j] = make_vec(x, y, z);
          }
        }

      }
    else if( _bReadFieldData && !sCurr.compare("FIELD") )
      {
      input >> sCurr;
      assert(!sCurr.compare("FieldData") );
      std::cout << "Found Field data...\n";
      input >> nCurr;
      _nNumOfFields = nCurr;
      assert(_nNumOfFields > 0);
      for( int nFieldCounter = 0; nFieldCounter < _nNumOfFields; ++nFieldCounter )
        {

        std::string name;
        input >> name;

        input >> nCurr;
        assert(nCurr == 1); // this version doesn't support multi dimensional fields

        input >> nCurr;
        assert(nCurr == _nNumOfPoints);

        input >> sCurr;
        assert(!sCurr.compare("float") );
        for( int i = 0; i < _nNumOfFibers; ++i )
          {
          int nLineLength = _lines[i].size();
          (*_fibers)[i].Fields[name].resize(nLineLength);
          for( int j = 0; j < nLineLength; ++j )
            {
            input >> (*_fibers)[i].Fields[name][j];
            }
          }
        }
      }
    }

  return 0; // should never reach eof
}
