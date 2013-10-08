/**
 * \file FiberFilter.cc
 * \brief Contains implementation of class FiberFilter
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#include "FiberFilter.h"
#include <iostream>
#include "fiber.h"

// MAIN FUNCTIONS

FiberFilter::FiberFilter()
{
  _bInSet = _bOutSet = _bRegionSet = false;
  _bConnectionMode = PASS;
  _bCalcInverse = false;
}

bool FiberFilter::Run()
{
  _lut = &_region->GetLut()[0];
  _refImage = Region::refImage;
  _nx = _region->SizeX();
  _ny = _region->SizeY();
  _nz = _region->SizeZ();
  _nSizeLUT = _region->GetLut().size();

  if( _bInSet && _bOutSet )
    {
    for( unsigned int i = 0; i < _inFibers->size(); i++ )
      {
      if( CheckConditions( (*_inFibers)[i]) )
        {
        _outFibers->push_back( (*_inFibers)[i]);
        }
      }
    return 0;
    }
  else
    {
    std::cout << "Error input or output fiber not set!\n";
    return 1;
    }
}

// SET FILTER OPTIONS //

void FiberFilter::SetInputFibers(const std::vector<Fiber> & fibers)
{
  _inFibers = &fibers;
  _bInSet = true;
}

void FiberFilter::SetOutputFibers(std::vector<Fiber> & fibers)
{
  _outFibers = &fibers;
  _bOutSet = true;
}

void FiberFilter::SetRegion(const Region & region)
{
  _region = &region;
  _bRegionSet = true;
}

void FiberFilter::SetConnectionMode(bool option)
{
  _bConnectionMode = option;
}

void FiberFilter::SetCalcInverse(bool option)
{
  _bCalcInverse = option;
}

// HELPER FUNCTIONS

bool FiberFilter::inRegion(const Fiber & fiber)
{

  typedef itk::Point<float, UIntImageType::ImageDimension> PointType;

  UIntImageType::IndexType index;
  PointType                point;
  int                      nMultIndex;

  std::vector<vec_t>::const_iterator cit;

  if( _bConnectionMode == PASS )
    {
    cit = fiber.Points.begin();
    }
  else if( _bConnectionMode == END )
    {
    cit = fiber.Points.end() - 3;
    }
  for( /* init above */; cit != fiber.Points.end(); ++cit )
    {

    point[0] = - (*cit)[0];
    point[1] = -(*cit)[1];
    point[2] = (*cit)[2];

    _refImage->TransformPhysicalPointToIndex(point, index);

    // clamping
    if( index[0] < 0 )
      {
      index[0] = 0;
      }
    if( index[1] < 0 )
      {
      index[1] = 0;
      }
    if( index[2] < 0 )
      {
      index[2] = 0;
      }
    if( index[0] > _nx - 1 )
      {
      index[0] = _nx - 1;
      }
    if( index[1] > _ny - 1 )
      {
      index[1] = _ny - 1;
      }
    if( index[2] > _nz - 1 )
      {
      index[2] = _nz - 1;
      }

    nMultIndex = index[0] * _ny * _nz + index[1] * _nz + index[2];
    for( unsigned int i = 0; i < _nSizeLUT; ++i )
      {
      if( nMultIndex == _lut[i] )
        {
        return 1;
        }
      }
    }
  return 0;
}

bool FiberFilter::CheckConditions(const Fiber & fiber)
{

  if( !_bCalcInverse )
    {
    return inRegion(fiber);
    }
  else if( _bCalcInverse )
    {
    return !inRegion(fiber);
    }

  std::cout << "Error: Should never reach this point!\n";
  exit(1);

}
