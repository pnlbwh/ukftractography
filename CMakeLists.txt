cmake_minimum_required(VERSION 3.5)

#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
include(ExternalProject)
include(ExternalProjectDependency)

#-----------------------------------------------------------------------------
if(APPLE)
  # Note: By setting CMAKE_OSX_* variables before any enable_language() or project() calls,
  #       we ensure that the bitness, and C++ standard library will be properly detected.
  include(SlicerBlockSetCMakeOSXVariables)
  mark_as_superbuild(
    VARS CMAKE_OSX_ARCHITECTURES:STRING CMAKE_OSX_SYSROOT:PATH CMAKE_OSX_DEPLOYMENT_TARGET:STRING
    ALL_PROJECTS
    )
endif()

set(PRIMARY_PROJECT_NAME UKFTractography)

## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  ${PRIMARY_PROJECT_NAME}_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with ${PRIMARY_PROJECT_NAME}_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the ${PRIMARY_PROJECT_NAME}, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
set(EXTENSION_NAME  UKFTractography)
set(EXTENSION_HOMEPAGE "http://www.nitrc.org/plugins/mwiki/index.php/ukftractography:MainPage")
set(EXTENSION_CATEGORY "Diffusion.Tractography")
set(EXTENSION_CONTRIBUTORS "Yogesh Rathi, Stefan Lienhard, Yinpeng Li, Martin Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner, Ryan Eckbo, Kent Williams, Hans Johnson, Peter Savadjiev, Carl-Fredrik Westin")
set(EXTENSION_MAJOR_VERSION 1)
set(EXTENSION_MINOR_VERSION 2)
set(EXTENSION_PATCH_VERSION 0)
set(EXTENSION_STATUS "Alpha")
set(EXTENSION_DEPENDS "NA")
set(EXTENSION_DESCRIPTION "This module traces fibers in a DWI Volume using the multiple tensor unscented Kalman Filter methology. ")
set(EXTENSION_LICENSE_SHORT_DESCRIPTION "Quadratic Programming Library QuadProg++ has LGPL or GPL")
set(EXTENSION_ICONURL "https://github.com/pnlbwh/ukftractography/raw/master/UKF_icon.png")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-CLIExtensionTemplate-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-CLIExtensionTemplate-screenshot-2.png")

#-----------------------------------------------------------------------------
option(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION "Build ${PRIMARY_PROJECT_NAME} as Slicer extension." OFF)

if (DEFINED Slicer_DIR)
  set(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION ON)
endif()

if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(default_superbuild_value OFF)
else()
  set(default_superbuild_value ON)
endif()

#-----------------------------------------------------------------------------
option(${PRIMARY_PROJECT_NAME}_SUPERBUILD "Build ${PRIMARY_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ${default_superbuild_value})
mark_as_advanced(${PRIMARY_PROJECT_NAME}_SUPERBUILD)

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${PRIMARY_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${PRIMARY_PROJECT_NAME})  # <- NOTE: Project name for pre-requisites is different form main project
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  # Building against Slicer
  project(${PRIMARY_PROJECT_NAME})             # <- NOTE: Here is the main project name setting
  include("${CMAKE_CURRENT_SOURCE_DIR}/${PRIMARY_PROJECT_NAME}.cmake")
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")

