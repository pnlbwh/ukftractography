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
#include "CompressedSensingCLP.h"
#include "CompressedSensing.h"
#include "DefaultGradients.h"
#include "nrrdIO.h"
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

/**
 * PrintMat
 * NOTE: during debugging, you can call this function to show contents
 * of Eigen matrices:
 * (gdb) call PrintMat(matrixName)
 */
#define DeclarPrintMat(type)                                    \
  void PrintMat(const type &mat)                                \
  {                                                             \
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");       \
    std::cerr << mat.format(CleanFmt) << std::endl;             \
  }

DeclarPrintMat(MatrixType)
DeclarPrintMat(Eigen::Vector3d)
DeclarPrintMat(Eigen::Vector3f)



bool GetGradientsFromFile(const std::string &filename, MatrixType &gradients)
{
  std::ifstream in(filename.c_str());
  // read whole file
  if(!in.is_open())
    {
    std::cerr << "Can't open " << filename << std::endl;
    return false;
    }

  in.seekg(0L, in.end);
  std::streampos length = in.tellg();
  in.seekg(0L, in.beg);
  char *buffer = new char[length];
  in.read(buffer, length);
  in.close();
  std::stringstream ss(buffer);

  bool commafile = std::strchr(buffer,',') != 0;

#if 0
  unsigned int pointcount = 0;

  std::vector<Eigen::Vector3d> pointlist;
  // if there is a comma in a file, assume csv file

  while(!ss.eof())
    {
    Eigen::Vector3d current;
    ss >> current(0);
    if(commafile)
      {
      ss.ignore(20,',');
      }
    ss >> current(1);
    if(commafile)
      {
      ss.ignore(20,',');
      }
    ss >> current(2);
    pointlist.push_back(current);
    pointcount++;
    }
  delete [] buffer;

  gradients.conservativeResize(3,pointcount);

  for(unsigned int i = 0; i < pointcount; ++i)
    {
    gradients.col(i) = pointlist[i];
    }
#else
  for(unsigned i = 0; !ss.eof(); ++i)
    {
    gradients.conservativeResize(3,i+1);
    Eigen::Vector3d current;
    ss >> current(0);
    if(commafile)
      {
      ss.ignore(20,',');
      }
    ss >> current(1);
    if(commafile)
      {
      ss.ignore(20,',');
      }
    ss >> current(2);
    gradients.col(i) = current;
    }
  delete [] buffer;
#endif
  return true;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  MatrixType newGradients;
  if(!replacementGradients.empty())
    {
    if(!GetGradientsFromFile(replacementGradients,newGradients))
      {
      return 1;
      }
    }
  else
    {
    SetDefaultGradients(newGradients);
    }
  return 0;
}
// function run_cs(input_dwi_fn,input_mask_fn,output_dwi_fn,new_gradients,saveC)
//
// fprintf('Starting testCS:\n    input_dwi_fn: %s\n    input_mask_fn: %s\n    output_dwi_fn: %s\n',input_dwi_fn,input_mask_fn,output_dwi_fn)
//
//
// if nargin >=4 && new_gradients
//     nu = new_gradients;
// else
//     nu = SetDefaultGradients();
// end
//
// [ rawDWI ] = nrrdLoadWithMetadata(input_dwi_fn);
// [ reformattedDWI ] = nrrdReformatAndNormalize(rawDWI);
// [ dwi_struct, metric, counts ] = BalanceDWIReplications( reformattedDWI );
// in_mask = nrrdLoadWithMetadata(input_mask_fn);
// mask = in_mask.data ~= 0;
//
// %save('/tmp/before_doCSestimate.mat','-v7.3')
// [estimatedSignal,estimatedGradients] = doCSestimate(dwi_struct, mask, nu);
//
// nrrdStrct = dwi_struct;
// nrrdStrct.data=single(estimatedSignal);
// nrrdStrct.gradientdirections=estimatedGradients;
//
// fprintf('Writing file to disk...\n');
//
// nrrdSaveWithMetadata(output_dwi_fn,nrrdStrct)
// if nargin == 5 && saveC
//   save(output_fn,'c','-v7.3');
// end
// end
//
// function [ nu ] =  SetDefaultGradients()
//     % icosahedron produces a antipodal symmetric set of points,
//     % so remove the redundant samples.
//
//     %ico1 = icosahedron(1); %42 gradient direction
//     ico2 = icosahedron(2); %162  <-- Good default for most cases
//     %ico3 = icosahedron(3); %642
//     nu = ico2; %42 gradient directions
//
//     n0 = size(nu,1)/2;
//     nu = nu(1:n0,:);
// end
