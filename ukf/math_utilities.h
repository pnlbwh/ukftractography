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
// Source: http://stackoverflow.com/questions/13690483/better-more-portable-method-of-defining-pi-in-c-c
//         http://www.geom.uiuc.edu/~huberty/math5337/groupe/digits.html

#  define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816406
#endif

#define UKF_PI M_PI

#define DEG_TO_RAD (UKF_PI/180.0)
#define RAD_TO_DEG (180.0/UKF_PI)
inline double DegToRad(const double deg) {
   return deg * DEG_TO_RAD;
}

inline double RadToDeg(const double rad) {
   return rad * RAD_TO_DEG;
}

#endif  // MATH_UTILITIES_H_
