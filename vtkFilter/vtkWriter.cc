/**
 * \file vtkWrite.cc
 * \brief Contains implementation of class vtkWriter
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#include "vtkWriter.h"
#include <iostream>
#include <fstream>

bool vtkWriter::Run()
{

  if( _fibers->empty() )
    {
    std::cout << "** No Fibers to write.";
    return 1;
    }

  std::ofstream vtkOutWriter;
  vtkOutWriter.open(_sOutputPath, std::ios::binary);

  if( !vtkOutWriter.is_open() )
    {
    std::cout << "Failed to open " << _sOutputPath << " for writing." << std::endl;
    return 1;
    }

  WriteHeader(vtkOutWriter);
  vtkOutWriter << "DATASET POLYDATA" << std::endl;
  WritePoints(vtkOutWriter);

  WriteLines(vtkOutWriter);
  vtkOutWriter << "POINT_DATA " << _nNumOfPoints << std::endl;
  // todo: tensors
  WriteFields(vtkOutWriter);

  vtkOutWriter.close();
  return 0;
}

void vtkWriter::SetInputFibers(std::vector<Fiber> & fibers)
{

  _fibers = &fibers;
  _nNumOfFibers = (*_fibers).size();

  _nNumOfFields = (*_fibers)[0].Fields.size();

  _fiberLengths.resize(_nNumOfFibers);
  _nNumOfPoints = 0;

  int nCurrSize;
  for( int i = 0; i < _nNumOfFibers; ++i )
    {
    nCurrSize = (*_fibers)[i].Points.size();
    _fiberLengths[i] = nCurrSize;
    _nNumOfPoints += nCurrSize;
    }
}

void vtkWriter::SetOutputPath(const std::string& path)
{

  _sOutputPath = path.c_str();

}

void vtkWriter::WriteHeader(std::ofstream & output)
{
  output << "#vtk DataFile Version 3.0" << std::endl;
  output << "Tracts filtered with vtkFilter" << std::endl;
  output << "ASCII" << std::endl;
}

void vtkWriter::WritePoints(std::ofstream & output)
{

  output << "POINTS " << _nNumOfPoints << " float";
  int nCounter = 0;
  for( int i = 0; i < _nNumOfFibers; ++i )
    {
    for( int j = 0; j < _fiberLengths[i]; ++j )
      {
      if( nCounter % 3 == 0 )
        {
        output << std::endl;
        }
      else
        {
        output << " ";
        }
      output << (*_fibers)[i].Points[j][0] << " "
             << (*_fibers)[i].Points[j][1] << " "
             << (*_fibers)[i].Points[j][2];
      nCounter++;
      }
    }
  output << std::endl;
}

void vtkWriter::WriteLines(std::ofstream & output)
{

  output << std::endl << "LINES ";
  output << _nNumOfFibers << " " << _nNumOfFibers + _nNumOfPoints << std::endl;

  int counter = 0;
  for( int i = 0; i < _nNumOfFibers; ++i )
    {
    output << _fiberLengths[i];
    for( int j = 0; j < _fiberLengths[i]; ++j )
      {
      output << " " << counter++;
      }
    output << std::endl;
    }
  output << std::endl;
}

void vtkWriter::WriteFields(std::ofstream & output)
{

  output << "FIELD FieldData " << _nNumOfFields << std::endl;
  Fiber::FieldMapType::const_iterator cit;
  int                                 nCounter;
  for( cit = (*_fibers)[0].Fields.begin(); cit != (*_fibers)[0].Fields.end(); ++cit )
    {
    nCounter = 0;

    output << cit->first;
    output << " 1 " << _nNumOfPoints;
    output << " float" << std::endl;
    for( int i = 0; i < _nNumOfFibers; ++i )
      {
      for( int j = 0; j < _fiberLengths[i]; ++j )
        {
        output << (*_fibers)[i].Fields[cit->first][j];
        nCounter++;
        if( nCounter % 9 == 0 )
          {
          output << std::endl;
          }
        else
          {
          output << " ";
          }
        }
      }

    output << std::endl;

    }

}
