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

/** Disable some common warnings in MS VC++ */
#if defined( _MSC_VER )

// 'identifier' : class 'type' needs to have dll-interface to be used by
// clients of class 'type2'
#pragma warning ( disable : 4251 )

#endif

#if defined(WIN32) && !defined(UKF_STATIC)
 #if defined(UKFBase_EXPORTS)
  #define UKFBASELIB_EXPORTS __declspec( dllexport )
 #else
  #define UKFBASELIB_EXPORTS __declspec( dllimport )
 #endif
#else
 #define UKFBASELIB_EXPORTS
#endif

#endif