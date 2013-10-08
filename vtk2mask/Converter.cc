/**
 * \file Converter.cc
 * \brief Contains implementations of Converter.h
*/

#include "Converter.h"
#include <iostream>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <fstream>
#include <cassert>

// CONSTRUCTOR //////////////////////////////////////////////////////////////////////////
Converter::Converter()
{
  // Changed to true if the respective variable is set
  _bOutputField        = false;
  _bHasLabel           = false;
  _bOutputVarianceMaps = false;
}

// RUN FUNCTION /////////////////////////////////////////////////////////////////////////
bool Converter::Run()
{

  if( _bVerbose )
    {
    std::cout << "** Reading DWI Data..\n";
    }
  CreateOutNrrd(_sReferenceVolume);

  if( _sLabelFile.empty() )
    {
    _bHasLabel = false;
    }
  else
    {
    _bHasLabel = true;
    LoadLabel();
    if( _bVerbose )
      {
      std::cout << "Using Label " << (uint)_nLabelOfInterest << " from label file\n";
      }
    }

  if( _bVerbose )
    {
    std::cout << "** Processing Fields...\n";
    }

  if( _bOutputField )
    {
    ProcessField(_sScalarFieldName, _sOutputVolume);
    }

  return 0;
}

// SETTERS ///////////////////////////////////////////////////////////////////////////////
void Converter::SetInputFibers(std::vector<Fiber> & fibers)
{
  _fibers = &fibers;
}

void Converter::SetReferenceFile(std::string & reference_volume)
{
  _sReferenceVolume = reference_volume;
}

void Converter::SetOutputVolumeFile(std::string & output_volume)
{
  _sOutputVolume = output_volume;
}

void Converter::SetStdDevFile(std::string & std_dev_volume)
{
  _sStandartDevVolume = std_dev_volume;
  _bOutputVarianceMaps = true;
}

void Converter::SetVerbose(bool verbose)
{
  _bVerbose = verbose;
}

void Converter::SetLabelFile(std::string & label_file)
{
  _sLabelFile = label_file;
}

void Converter::SetLabelNumber(int label)
{
  _nLabelOfInterest = label;
  _bHasLabel = true;
}

void Converter::SetFieldName(std::string & field_name)
{
  _sScalarFieldName =  field_name;
  _bOutputField = true;
}

// PRIVATE MEMBER FUNCTIONS ////////////////////////////////////////////////////////////////
// Distribute the Work
void Converter::ProcessField(const std::string & field_name, const std::string & path)
{
  _matField = new Array3D<std::vector<float> >(_nDimX, _nDimY, _nDimZ);
  _nrrdDataOut->FillBuffer( 0 ); // start a new scalar map
  if( _bOutputVarianceMaps )
    {
    _nrrdVarDataOut->FillBuffer( 0 );
    }
  FillMatField(field_name);
  AverageVoxels();
  WriteOutNrrd(path);
  delete _matField;
}

// Fill the matrix of vectors with the values
bool Converter::FillMatField(const std::string & field_name)
{

  int            x, y, z;
  int            ec = 0;
  FloatPointType point;

  for( unsigned int i = 0; i < _fibers->size(); ++i )
    {
    int nNumOfPoints = (*_fibers)[i].Points.size();
    for( int j = 0; j < nNumOfPoints; ++j )
      {

      // ITK uses a different measurement frame than VTK. The flipping of x and y
      // might have something todo with that. Not sure yet, though.
      point[0] = -(*_fibers)[i].Points[j][0];
      point[1] = -(*_fibers)[i].Points[j][1];
      point[2] =  (*_fibers)[i].Points[j][2];

      FloatImageType::IndexType pixelIndex;
      _nrrdDataOut->TransformPhysicalPointToIndex(point, pixelIndex);

      x = pixelIndex[0];
      y = pixelIndex[1];
      z = pixelIndex[2];

      // Take care of rounding errors at borders of volume
      if( x < 0 )
        {
        ec++;
        continue;
        }
      if( y < 0 )
        {
        ec++;
        continue;
        }
      if( z < 0 )
        {
        ec++;
        continue;
        }
      if( x >= _nDimX )
        {
        ec++;
        continue;
        }
      if( y >= _nDimY )
        {
        ec++;
        continue;
        }
      if( z >= _nDimZ )
        {
        ec++;
        continue;
        }
      if( !field_name.empty() )
        {
        _matField->_[x][y][z].push_back( (*_fibers)[i].Fields[field_name][j]);
        }
      else
        {
        _matField->_[x][y][z].push_back(1); // fiber passed through this point
        }
      }
    }
  if( _bVerbose )
    {
    if( ec > 0 )
      {
      std::cout << "Number of clipped points: " << ec << std::endl;
      }
    }
  return 0;
}

// Average the vectors in the matrix
void Converter::AverageVoxels()
{
  FloatImageType::IndexType current;
  float                     mean, stddev;
  int                       counter = 0;

  for( int x = 0; x < _nDimX; ++x )
    {
    for( int y = 0; y < _nDimY; ++y )
      {
      for( int z = 0; z < _nDimZ; ++z )
        {
        current[0] = x;
        current[1] = y;
        current[2] = z;
        if( (!_bHasLabel || _Label->GetPixel( current) == (uint)_nLabelOfInterest) && _matField->_[x][y][z].size() > 0 )
          {
          mean = CalcMean(_matField->_[x][y][z]);
          _nrrdDataOut->SetPixel( current, mean );
          if( _bOutputVarianceMaps )
            {
            stddev = CalcStdDev( _matField->_[x][y][z], mean );
            _nrrdVarDataOut->SetPixel( current, stddev);
            }
          counter++;
          }
        }
      }
    }
}

// HELPER FUNCTIONS ///////////////////////////////////////////////////////////////
// file I/O
void Converter::LoadLabel()
{
  std::cout << "** Reading Labels...\n";

  typedef itk::ImageFileReader<UIntImageType> FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(_sLabelFile);
  reader->Update();
  _Label = reader->GetOutput();
}

// Create the output Volume
void Converter::CreateOutNrrd(const std::string & nrrd_path)
{

  typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;

  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(nrrd_path);
  reader->Update();
  DiffusionImageType::Pointer refNrrd = reader->GetOutput();

  // Copy information to outfile
  FloatImageType::IndexType start;
  start[0] = start[1] = start[2] = 0;

  FloatImageType::SizeType size = refNrrd->GetLargestPossibleRegion().GetSize();

  _nDimX = size[0];
  _nDimY = size[1];
  _nDimZ = size[2];

  FloatImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  _nrrdDataOut = FloatImageType::New();
  _nrrdDataOut->SetSpacing( refNrrd->GetSpacing() );
  _nrrdDataOut->SetOrigin( refNrrd->GetOrigin() );
  _nrrdDataOut->SetDirection( refNrrd->GetDirection() );
  _nrrdDataOut->SetRegions( region );
  _nrrdDataOut->Allocate();

  if( _bOutputVarianceMaps )
    {
    _nrrdVarDataOut = FloatImageType::New();
    _nrrdVarDataOut->SetSpacing( refNrrd->GetSpacing() );
    _nrrdVarDataOut->SetOrigin( refNrrd->GetOrigin() );
    _nrrdVarDataOut->SetDirection( refNrrd->GetDirection() );
    _nrrdVarDataOut->SetRegions( region );
    _nrrdVarDataOut->Allocate();
    }
}

// Write the output Volume
void Converter::WriteOutNrrd(const std::string & out_path)
{

  typedef itk::ImageFileWriter<FloatImageType> WriterType;

  itk::NrrdImageIO::Pointer io = itk::NrrdImageIO::New();
  // io->SetFileType( itk::ImageIOBase::ASCII );
  WriterType::Pointer nrrdWriter = WriterType::New();

  nrrdWriter->SetInput( _nrrdDataOut );
  nrrdWriter->SetFileName(out_path);

  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject e )
    {
    std::cout << e << std::endl;
    }

  if( _bOutputVarianceMaps )
    {

    std::stringstream sVarFilePath;

    nrrdWriter->SetInput( _nrrdVarDataOut );
    nrrdWriter->SetFileName(_sStandartDevVolume);

    try
      {
      nrrdWriter->Update();
      }
    catch( itk::ExceptionObject e )
      {
      std::cout << e << std::endl;
      }
    }

}

// Tools
float Converter::CalcStdDev(const std::vector<float> & vec, const float & mean)
{
  float sum = 0;

  for( uint i = 0; i < vec.size(); ++i )
    {
    sum += (vec[i] - mean) * (vec[i] - mean);
    }
  return sqrt(sum / vec.size() );
}

float Converter::CalcMean(const std::vector<float> &  vec)
{
  return std::accumulate( vec.begin(), vec.end(), 0.0f ) / vec.size();
}
