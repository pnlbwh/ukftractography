include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

#-----------------------------------------------------------------------------
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#-----------------------------------------------------------------------------
find_package(ITK REQUIRED)
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # Incorporate with Slicer nicely
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
if(NOT Eigen_INCLUDE_DIR)
  if(NOT Eigen_DIR)
    message(FATAL_ERROR "Missing Eigen_DIR path, can't find Eigen library includes")
  endif()
  set(Eigen_INCLUDE_DIR
    ${Eigen_DIR}/../Eigen)
endif()
include_directories(${Eigen_INCLUDE_DIR})

find_package(ZLIB REQUIRED)

#
#-----------------------------------------------------------------------------
find_package(Teem REQUIRED)
include(${Teem_USE_FILE})
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(TEEM_LIB teem)
else()
  #
  # due to forcing all dependcy builds to put outputs into lib and bin at the
  # top level build directory, the Teem_LIBRARY_DIRS var is wrong; have to add
  # top level build dir here
  find_library(TEEM_LIB teem PATHS ${CMAKE_CURRENT_BINARY_DIR}/../lib)
  message("TEEM_LIB:${TEEM_LIB}")
endif()

#-----------------------------------------------------------------------------
add_subdirectory(common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
add_subdirectory(ukf)

if(NOT ${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  option(USE_fibertractdispersion "Build the fibertractdispersion program" OFF)
  if(USE_fibertractdispersion)
    add_subdirectory(fibertractdispersion)
  endif()
  option(USE_CompressedSensing "Build the CompressedSensing program" OFF)
  if(USE_CompressedSensing)
    add_subdirectory(CompressedSensing)
  endif()
  add_subdirectory(vtk2mask)
  add_subdirectory(vtkFilter)
endif()

#-----------------------------------------------------------------------------
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  include(${Slicer_EXTENSION_CPACK})
endif()

