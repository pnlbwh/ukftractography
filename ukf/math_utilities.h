/**
* \file math_utilities.h
* \brief Implements some math functions needed in various parts of the code
*
* This file reimplements some functionality from math.h in order to overcome some cross-platform problems.
* There are no execution speed repercussions resulting from this substitution.
*
* \author Christian Baumgartner (c.f.baumgartner@gmail.com)
*/

#ifndef MATH_UTILITIES_H_
#define MATH_UTILITIES_H_

#include <limits>
#ifndef M_PI
/** Approximation of Pi, if necesseary */
#define M_PI 3.14159265
#endif

/**
* \brief Returns true for a variable of every type that has infinity defined in 
* std and is larger than its maximal value. 
* \param value A variable of any of the template types.
*/
template<typename T>
inline bool isinf(T value)
{
  return std::numeric_limits<T>::has_infinity && value == std::numeric_limits<T>::infinity();
}

/**
* \brief Returns true if the given value is not a number (NaN).
* \param value A value of any of the template types.
*/
template<typename T>
inline bool isnan(T value)
{
  return value != value;
}

#endif  // MATH_UTILITIES_H_
