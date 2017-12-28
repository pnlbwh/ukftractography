/*=auto=========================================================================

 Portions (c) Copyright 2005 Brigham and Women's Hospital (BWH)
 All Rights Reserved.

 See COPYRIGHT.txt
 or http://www.slicer.org/copyright/copyright.txt for details.

 Program:   3D Slicer

=========================================================================auto=*/


// .NAME __UKFExport - manage Windows system differences
// .SECTION Description
// The __UKFExport captures some system differences between Unix
// and Windows operating systems.

#ifndef __UKFExport_h
#define __UKFExport_h

#if defined(WIN32) && !defined(UKF_STATIC)
 #if defined(UKF_EXPORTS)
  #define UKFTRACTOGRAPHYLIB_EXPORT __declspec( dllexport )
 #else
  #define UKFTRACTOGRAPHYLIB_EXPORT __declspec( dllimport )
 #endif
#else
 #define UKFTRACTOGRAPHYLIB_EXPORT
#endif
#endif
