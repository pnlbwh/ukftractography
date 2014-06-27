include(${CMAKE_CURRENT_LIST_DIR}/Common.cmake)


#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
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
if(Slicer_BUILD_${PROJECT_NAME})
  set(ITK_NO_IO_FACTORY_REGISTER_MANAGER 1) # Incorporate with Slicer nicely
endif()
include(${ITK_USE_FILE})

if(NOT Eigen_INCLUDE_DIR)
  if(NOT Eigen_DIR)
    message(FATAL_ERROR "Missing Eigen_DIR path, can't find Eigen library includes")
  endif()
  set(Eigen_INCLUDE_DIR
    ${Eigen_DIR}/../Eigen)
endif()

if(NOT Teem_FOUND)
find_package(Teem REQUIRED)
include(${Teem_USE_FILE})
endif()

find_library(TEEM_LIB name teem PATHS ${Teem_LIBRARY_DIRS})
if("${TEEM_LIB}" EQUAL "TEEM_LIB-NOTFOUND")
  message(FATAL_ERROR "Can't find Teem library TEEM_LIB")
else()
  message("TEEM_LIB=${TEEM_LIB}")
endif()


include_directories(${Eigen_INCLUDE_DIR})
#-----------------------------------------------------------------------------
add_subdirectory(common)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/common)
add_subdirectory(ukf)
add_subdirectory(fibertractdispersion)
add_subdirectory(CompressedSensing)
add_subdirectory(vtk2mask)
add_subdirectory(vtkFilter)
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()

