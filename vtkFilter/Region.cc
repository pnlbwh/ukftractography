/**
 * \file Region.cc
 * \brief Contains implementations of class Region
 * \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#include "Region.h"
#include <iostream>

Region::UIntImageType::Pointer Region::refImage;

Region::Region(const std::string & path, const std::vector<int> & labels)
{

  if (path.empty()) {
    _bEmpty = true;
  } else {
    _bEmpty = false;

    typedef itk::ImageFileReader<UIntImageType> FileReaderType;
    typedef itk::ImageRegionConstIterator<UIntImageType> ConstIteratorType;

    FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName(path);
    reader->Update();
    refImage = reader->GetOutput();

    UIntImageType::RegionType region;
    UIntImageType::IndexType index;
    UIntImageType::SizeType size;

    region = refImage->GetLargestPossibleRegion();
    size = region.GetSize();

    _nSizeX = size[0];
    _nSizeY = size[1];
    _nSizeZ = size[2];

    ConstIteratorType cit( refImage, region );

    while( !cit.IsAtEnd() ) {
      if (inVector(cit.Get(), labels)) {
        index = cit.GetIndex();
        _lut.push_back(_nSizeZ * _nSizeY * index[0] + _nSizeZ * index[1] + index[2]);
      }
      ++cit;
    }

  }

}

bool Region::inVector(const int & n, const std::vector<int> & v)
{
  for (unsigned int i = 0; i < v.size(); ++i) {
    if (v[i] == n)
      return true;
  }
  return false;
}