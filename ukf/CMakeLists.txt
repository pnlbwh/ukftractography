cmake_minimum_required(VERSION 2.8.7)


#-----------------------------------------------------------------------------
set(MODULE_NAME UKFTractography)

#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR AND NOT UKFTractography_SuperBuild)
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
endif()

if (UKFTractography_SuperBuild)
    find_package(SlicerExecutionModel REQUIRED)
    include(${SlicerExecutionModel_USE_FILE})
    if (NOT Teem_FOUND)
        find_package(Teem REQUIRED)
    endif (NOT Teem_FOUND)
    include(${Teem_USE_FILE})

    if (NOT ITK_FOUND)
        find_package(ITK REQUIRED)
    endif (NOT ITK_FOUND)
    include(${ITK_USE_FILE})

    if (NOT GenerateCLP_FOUND)
      find_package(GenerateCLP REQUIRED)
    endif (NOT GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})

    #if (NOT Boost_FOUND)
      #find_package(Boost 1.41.0 COMPONENTS thread REQUIRED)
    #endif (NOT Boost_FOUND)

    link_libraries(${Boost_LIBRARY_DIRS}/libboost_thread.a)
endif()

#-----------------------------------------------------------------------------
add_subdirectory(Testing)

#-----------------------------------------------------------------------------


set(MODULE_INCLUDE_DIRECTORIES
    ${Boost_INCLUDE_DIRS}
    #${Teem_INCLUDE_DIRS}
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
    #${Boost_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
  )


# In the config file we integrate a boost toggle
configure_file (
 "${PROJECT_SOURCE_DIR}/config.h.in"
 "${PROJECT_SOURCE_DIR}/config.h"
 )


if (NOT GenerateCLP_FOUND)
  find_package(GenerateCLP REQUIRED)
endif (NOT GenerateCLP_FOUND)
include(${GenerateCLP_USE_FILE})
#-----------------------------------------------------------------------------
message(${Boost_LIBRARY_DIRS})
if(UKFTractography_SuperBuild)
    SEMMacroBuildCLI(
        NAME ${MODULE_NAME}
        LINK_DIRECTORIES ${Boost_LIBRARY_DIRS}
        TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
        INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
        ADDITIONAL_SRCS ${MODULE_SRCS}
        EXECUTABLE_ONLY
        )
else()
    SEMMacroBuildCLI(
        NAME ${MODULE_NAME}
        TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
        INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
        ADDITIONAL_SRCS ${MODULE_SRCS}
        )
endif()
