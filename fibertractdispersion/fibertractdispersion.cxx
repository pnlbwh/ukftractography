/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile$
Language:  C++
Date:      $Date: 2007-08-31 11:20:20 -0500 (Fri, 31 Aug 2007) $
Version:   $Revision: 10358 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "fibertractdispersionCLP.h"
#include <sstream>
#include <map>
#include <vector>
#include "fiber.h"
#include "fiberbundle.h"
#include "computedispersion.h"

int main( int argc, char * argv[] )
{
  // Default command line  ./fibertractdispersion CompressedSensingDWI/TestSuite/disp_scale1_numdirs10_1_1_circles.vtk test.vtk
  // Ensure that the test.vtk file can be opened in Slicer, and colored by "TEST_ATTRIBUTE".

  PARSE_ARGS;

  /*  PART 1:  Read in the input fibers */
  fiberbundle  mybundle;
  mybundle.ReadFibers( inputFiberBundle );

  /* PART 2: Algorithm stuff here */
  computedispersion(mybundle,dispersionScale,numberOfSamplingDirections,"",tractSubsampling,fiberPointSubsampling);

  /* PART 3: Write bundle to disk */
  mybundle.WriteFibers( outputFiberBundle,writeAscii,writeUnCompressed );

  return 0;
}
