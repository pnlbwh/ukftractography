
cmake_minimum_required(VERSION 2.8.7)


#-----------------------------------------------------------------------------
set(MODULE_NAME fibertractdispersion)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
    ${CMAKE_CURRENT_SOURCE_DIR}/../ukf
 )

include_directories( ${MODULE_INCLUDE_DIRECTORIES} )

set(MODULE_SRCS
  fibertractdispersion.cxx
  computedispersion.cxx
  fiberbundle.cxx
  )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${VTK_LIBRARIES}
    ${TEEM_LIB}
  )

#-----------------------------------------------------------------------------
enable_testing()
add_subdirectory(Testing)

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )
