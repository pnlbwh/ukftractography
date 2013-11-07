/**
 * \file fiber.h
 * \brief Contains a class for fiber representation
*/

#ifndef FIBER_H_
#define FIBER_H_

#include <vector>
#include "linalg.h"
#include <iostream>
#include <map>

/**
 * \struct Fiber
 * \brief Describes a single fiber
*/
class Fiber
  {

public:

  /** A map for accessing tensors. Key the tensor name; value a vector of 3D matrices. */
  typedef std::map<std::string, stdMat_t, std::less<std::string>,
         Eigen::aligned_allocator<std::pair<const std::string, stdMat_t> > > TensorMapType;


  /** A map with the fieldname as key and a 2D ukfPrecisionType vector for the field contents */
  typedef std::map<std::string, std::vector<float> > FieldMapType;

  /** The points in RAS space of the fiber */
  stdVec_t Points;

  /** The fields of the fiber */
  FieldMapType Fields;

  /** The tensors of the fiber */
  TensorMapType Tensors;

  };

#endif
