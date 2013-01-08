# ---------- Build script for Unscented Kalman Filter Tractography ---------------------------
# Additional Information:
# -- Homepage: http://www.nitrc.org/plugins/mwiki/index.php/ukftractography:MainPage
# -- Contributors: Yogesh Rathi, Stefan Lienhard, Yinpeng Li, Martin Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner
# -- Version: 1.2
# -- Status: alpha
# --------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
# Settings and Options
#-----------------------------------------------------------------------------

include(BuildTools/BuildMacros.cmake)
include(BuildTools/PreventInSourceBuilds.cmake)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

string(TOUPPER CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} BUILD_C_FLAGS)
string(TOUPPER CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} BUILD_CXX_FLAGS)

set(BUILD_SHARED_LIBS OFF)
set(EXECUTABLE_OUTPUT_PATH bin)

option(BUILD_TESTING "Build testing for this project" ON)
option(BUILD_OPTIONAL_TOOLS "Build the two optional tools vtkFilter, and vtk2mask for this project" ON)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

enable_testing()

#-----------------------------------------------------------------------------
# Module description
#-----------------------------------------------------------------------------
set(EXTENSION_NAME "UKFTractography")

#-----------------------------------------------------------------------------
# Find Packages
#-----------------------------------------------------------------------------

if (NOT UKFTractography_SuperBuild)
  
  GETSLICERVERSION(Slicer_VERSION)
  message(STATUS "Slicer Version: ${Slicer_VERSION}") # 0 means slicer not found

endif (NOT UKFTractography_SuperBuild)

# A Warning for people trying to build without Slicer and without SuperBuild
if (Slicer_VERSION EQUAL 0 AND NOT UKFTractography_SuperBuild)
  message("----------------------------------------------------------------------------------------------------")
  message("If you are having trouble to compile the Module without Slicer this way, use the SuperBuild instead.")
  message("This is, run the following commands from your build directory:")
  message("   cmake <path-to-ukfsource>/SuperBuild")
  message("   make")
  message("It will take care of all dependancies for you")
  message("If you have Slicer, set the Slicer3_DIR (Sl3) or the Slicer_DIR (Sl4) when runing cmake. Example:")
  message("   cmake -DSlicer_DIR:PATH=<path-to-your-slicer>/Slicer-build/ <path-to-ukfsource>")
  message("   make")
  message("----------------------------------------------------------------------------------------------------")
endif (Slicer_VERSION EQUAL 0 AND NOT UKFTractography_SuperBuild)

if (NOT ITK_FOUND)
  find_package(ITK REQUIRED)
endif (NOT ITK_FOUND)
include(${ITK_USE_FILE})

if (NOT Teem_FOUND)
  find_package(Teem REQUIRED)
endif (NOT Teem_FOUND)
include(${Teem_USE_FILE})

if (NOT GenerateCLP_FOUND)
  find_package(GenerateCLP REQUIRED)
endif (NOT GenerateCLP_FOUND)
include(${GenerateCLP_USE_FILE})

if (NOT Boost_FOUND)
  find_package(Boost 1.41.0 COMPONENTS thread REQUIRED)
endif (NOT Boost_FOUND)

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

#-----------------------------------------------------------------------------
# Set sources etc.
#-----------------------------------------------------------------------------

# Additional includes
set(cli_module_include_directories
    ${Boost_INCLUDE_DIRS}
  )

# Source files
set(cli_module_SRCS 
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
  
# Additional Target libraries
set(cli_module_target_libraries
    ${ITK_LIBRARIES}
    teem 
    ${Boost_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
  )

# In the config file we integrate a boost toggle
configure_file (
 "${PROJECT_SOURCE_DIR}/config.h.in"
 "${PROJECT_SOURCE_DIR}/config.h"
 )

#-----------------------------------------------------------------------------
# Building
#-----------------------------------------------------------------------------

GENERATECLPMACRO(
  MAIN_SRC "${EXTENSION_NAME}.cxx"
  ADDITIONAL_SRCS ${cli_module_SRCS}
  XML_FILE "${EXTENSION_NAME}.xml"
  EXECUTABLE_NAME ${EXTENSION_NAME}
  TARGET_LIBRARIES ${cli_module_target_libraries}
)   

#-----------------------------------------------------------------------------
# Build additional tools 
#-----------------------------------------------------------------------------
if (BUILD_OPTIONAL_TOOLS)
  add_subdirectory(OptionalTools)
endif (BUILD_OPTIONAL_TOOLS)

#-----------------------------------------------------------------------------
# Installing
#-----------------------------------------------------------------------------
if(UKFTractography_SuperBuild)
  install(TARGETS ${EXTENSION_NAME} RUNTIME DESTINATION bin)
  # Optional tools. There were problems using install TARGETS from respective sub cmake files.
  install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/vtk2mask DESTINATION bin) 
  install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/vtkFilter DESTINATION bin) 
endif(UKFTractography_SuperBuild)

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

