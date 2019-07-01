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
find_package(ITK REQUIRED)
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # Incorporate with Slicer nicely
endif()
include(${ITK_USE_FILE})

#-----------------------------------------------------------------------------
set(VTK_FOUND OFF)
find_package(VTK COMPONENTS
      vtkCommonSystem
      vtkCommonCore
      vtkCommonSystem
      vtkCommonMath
      vtkCommonMisc
      vtkCommonTransforms
      vtkIOLegacy
      vtkIOXML
      REQUIRED)
if(VTK_USE_FILE)
  include(${VTK_USE_FILE})
endif()

#-----------------------------------------------------------------------------
if(DEFINED Eigen_INCLUDE_DIR)
  include_directories(${Eigen_INCLUDE_DIR})
else()
  if(DEFINED Eigen_DIR)
    set(Eigen_INCLUDE_DIR
      ${Eigen_DIR}/../Eigen)
    include_directories(${Eigen_INCLUDE_DIR})
  else()
    set (Eigen_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen)
    set (Eigen_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/Eigen-build)
    ExternalProject_Add(
      UKF-Eigen
      DOWNLOAD_DIR      ${Eigen_DIR}
      SOURCE_DIR        ${Eigen_DIR}
      BINARY_DIR        ${Eigen_BUILD_DIR}
      GIT_REPOSITORY    ${Eigen_GIT_REPOSITORY}
      GIT_TAG           ${Eigen_GIT_TAG}
      CONFIGURE_COMMAND ""
      BUILD_COMMAND     ""
      INSTALL_COMMAND   ""
      )
    include_directories(${Eigen_DIR})
  endif()
endif()

find_package(ZLIB REQUIRED)

#
#-----------------------------------------------------------------------------
find_package(Teem REQUIRED)
include(${Teem_USE_FILE})
if(NOT ${PRIMARY_PROJECT_NAME}_SUPERBUILD)
  set(TEEM_LIB teem)
else()
  # due to forcing all dependency builds to put outputs into lib and bin at the
  # top level build directory, the Teem_LIBRARY_DIRS var is wrong; have to add
  # top level build dir here
  find_library(TEEM_LIB teem PATHS ${CMAKE_CURRENT_BINARY_DIR}/../lib)
  message("TEEM_LIB:${TEEM_LIB}")
endif()

#-----------------------------------------------------------------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
set(UKF_STATIC)
if(NOT BUILD_SHARED_LIBS)
    set(UKF_STATIC 1)
    if(WIN32)
        add_definitions("-DUKF_STATIC")
    endif()
endif()
add_subdirectory(ukf)
add_subdirectory(UKFTractography)

if(NOT ${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  option(USE_fibertractdispersion "Build the fibertractdispersion program" ON)
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

if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  add_subdirectory(InteractiveUKF)
endif()

#-----------------------------------------------------------------------------
if(${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  include(${Slicer_EXTENSION_GENERATE_CONFIG})
  include(${Slicer_EXTENSION_CPACK})
endif()

