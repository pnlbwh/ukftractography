cmake_minimum_required(VERSION 2.8.7)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  set(EXTENSION_NAME  UKFTractography)
  set(EXTENSION_HOMEPAGE "http://www.nitrc.org/plugins/mwiki/index.php/ukftractography:MainPage")
  set(EXTENSION_CATEGORY "Diffusion.Tractography")
  set(EXTENSION_CONTRIBUTORS "Yogesh Rathi, Stefan Lienhard, Yinpeng Li, Martin Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner")
  set(EXTENSION_MAJOR_VERSION 1)
  set(EXTENSION_MINOR_VERSION 2)
  set(EXTENSION_PATCH_VERSION 0)
  set(EXTENSION_STATUS "Alpha")
  set(EXTENSION_DESCRIPTION "This module traces fibers in a DWI Volume using the multiple tensor unscented Kalman Filter methology. ")
  set(EXTENSION_LICENSE_SHORT_DESCRIPTION "Quadratic Programming Library QuadProg++ has LGPL or GPL")
  set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/CLIExtensionTemplate/CLIExtensionTemplate.png?revision=19437&view=co")
  set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-CLIExtensionTemplate-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-CLIExtensionTemplate-screenshot-2.png")
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
else()
  #-----------------------------------------------------------------------------
  find_package(SlicerExecutionModel REQUIRED GenerateCLP)
  include(${GenerateCLP_USE_FILE})
  include(${SlicerExecutionModel_USE_FILE})
  include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)
endif()

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
#-----------------------------------------------------------------------------
add_subdirectory(common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
add_subdirectory(ukf)
add_subdirectory(vtk2mask)
add_subdirectory(vtkFilter)
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

