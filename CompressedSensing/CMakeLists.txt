
cmake_minimum_required(VERSION 2.8.7)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
message("Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
#-----------------------------------------------------------------------------
set(MODULE_NAME CompressedSensing)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR AND NOT UKFTractography_SuperBuild)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

if(NOT SlicerExecutionModel_FOUND)
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
endif()

if(NOT Teem_FOUND)
find_package(Teem REQUIRED)
include(${Teem_USE_FILE})
endif()

if(NOT ITK_FOUND)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
endif()

if(NOT GenerateCLP_FOUND)
find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})
endif()

if(NOT VTK_FOUND)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
endif()

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
    ${Teem_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
 )

include_directories( ${MODULE_INCLUDE_DIRECTORIES} )

set(MODULE_SRCS
  CompressedSensing.cxx
  SphericalHarmonics.cxx
  DefaultGradients.cxx
  DoCSEstimate.cxx
  icosahedron3.cxx
  )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    teem
  )

#-----------------------------------------------------------------------------
enable_testing()
add_subdirectory(Testing)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES} ${VTK_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )
