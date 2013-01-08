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
set(MODULE_NAME UKFTractography)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
#add_subdirectory(Logic)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
    ${Boost_INCLUDE_DIRS}
  )

set(MODULE_SRCS
    tractography.cc 
    unscented_kalman_filter.cc 
    filter_model.cc   
    seed.cc 
    fiber.cc 
    NrrdData.cc 
    vtk_writer.cc 
    dwi_normalize.cc 
    timer.cc
    utilities.cc 
    thread.cc 
    QuadProg++_vnl.cc
    )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    teem 
    ${Boost_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

