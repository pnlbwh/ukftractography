cmake_minimum_required(VERSION 2.8.7)


#-----------------------------------------------------------------------------
set(MODULE_NAME UKFTractography)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

#-----------------------------------------------------------------------------
add_subdirectory(Testing)

#-----------------------------------------------------------------------------
if (NOT Boost_FOUND)
  find_package(Boost 1.41.0 COMPONENTS thread REQUIRED)
endif (NOT Boost_FOUND)

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
