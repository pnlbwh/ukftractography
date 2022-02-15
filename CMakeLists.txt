# ==== Define cmake build policies that affect compilation and linkage default behaviors
#
# Set the NEWEST_VALIDATED_POLICIES_VERSION string to the newest cmake version
# policies that provide successful builds. By setting NEWEST_VALIDATED_POLICIES_VERSION
# to a value greater than the oldest policies, all policies between
# OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION (used for this build)
# are set to their NEW behaivor, thereby suppressing policy warnings related to policies
# between the OLDEST_VALIDATED_POLICIES_VERSION and CMAKE_VERSION.
#
# CMake versions greater than the NEWEST_VALIDATED_POLICIES_VERSION policies will
# continue to generate policy warnings "CMake Warning (dev)...Policy CMP0XXX is not set:"
#
set(OLDEST_VALIDATED_POLICIES_VERSION "3.16.3")
set(NEWEST_VALIDATED_POLICIES_VERSION "3.19.7")
cmake_minimum_required(VERSION ${OLDEST_VALIDATED_POLICIES_VERSION}...${NEWEST_VALIDATED_POLICIES_VERSION} FATAL_ERROR)

#-----------------------------------------------------------------------------
# Enable C++14
#-----------------------------------------------------------------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "${_msg} - C++${CMAKE_CXX_STANDARD}")
if(NOT CMAKE_CXX_STANDARD MATCHES "^(14|17|20)$")
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be set to 14, 17, or 20")
endif()


#-----------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/CMake)
include(SlicerBlockSetCMakeOSXVariables)


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
# Extension meta-information
#-----------------------------------------------------------------------------
set(EXTENSION_NAME  UKFTractography)
set(EXTENSION_HOMEPAGE "http://www.nitrc.org/plugins/mwiki/index.php/ukftractography:MainPage")
set(EXTENSION_CATEGORY "Diffusion.Tractography")
set(EXTENSION_CONTRIBUTORS "Yogesh Rathi, Stefan Lienhard, Yinpeng Li, Martin Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner, Ryan Eckbo, Kent Williams, Hans Johnson, Peter Savadjiev, Carl-Fredrik Westin")
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
# Additional CXX/C Flags
#-----------------------------------------------------------------------------
set(ADDITIONAL_C_FLAGS "" CACHE STRING "Additional C Flags")
mark_as_advanced(ADDITIONAL_C_FLAGS)
set(ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional CXX Flags")
mark_as_advanced(ADDITIONAL_CXX_FLAGS)

#-----------------------------------------------------------------------------
# Required by Slicer extension build system for reasons. See
#    https://github.com/Slicer/Slicer/commit/b160ec13f276a86306513954ef8b08a5332afc2e
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};RuntimeLibraries;/")

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------

if(${PRIMARY_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${PRIMARY_PROJECT_NAME} LANGUAGES C CXX)  # <- NOTE: Project name for pre-requisites is different form main project
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  # Building against Slicer
  project(${PRIMARY_PROJECT_NAME} LANGUAGES C CXX)             # <- NOTE: Here is the main project name setting
  include("${CMAKE_CURRENT_SOURCE_DIR}/${PRIMARY_PROJECT_NAME}.cmake")
  return()
endif()
message(FATAL_ERROR "You should never reach this point !")
