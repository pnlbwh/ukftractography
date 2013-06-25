/**
 * \file Array3D.h
 * \brief Contains class for large 3D Arrays
*/

#ifndef ARRAY3D_H_
#define ARRAY3D_H_

/**
 * \class Array3D
 * \brief Template Class for large 3D Arrays
 *
 * Note: The vector of vector of vector notation for 3D arrays has problems with such large dimensions
*/
template <class T>
class Array3D
{
private:

  /** X Dimension */
  int _nDimX;

  /** Y Dimension */
  int _nDimY;

  /** Z Dimension */
  int _nDimZ;
public:

  /** Pointer to the data */
  T * * *_;

  /** Constructor, allocates the space in each dimension */
  Array3D(const int x, const int y, const int z)
  {
    _nDimX = x;
    _nDimY = y;
    _nDimZ = z;

    _ = new T * *[_nDimX];
    for( int i = 0; i < _nDimX; ++i )
      {
      _[i] = new T *[_nDimY];
      for( int j = 0; j < _nDimY; ++j )
        {
        _[i][j] = new T[_nDimZ];
        }
      }
  }

  /** Deconstructor, deletes the data */
  ~Array3D()
  {
    delete _; //HACK: This is a huge memory leak here
  };                          // maybe have to iterate through matrix and delete

  /** Get the X dimensions */
  int dimX()
  {
    return _nDimX;
  }

  /** Get the Y dimensions */
  int dimY()
  {
    return _nDimY;
  }

  /** Get the Z dimensions */
  int dimZ()
  {
    return _nDimZ;
  }

};

#endif
