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

    find_package(Teem REQUIRED)
    include(${Teem_USE_FILE})

    find_package(ITK REQUIRED)
    include(${ITK_USE_FILE})

    find_package(GenerateCLP REQUIRED)
    include(${GenerateCLP_USE_FILE})

endif()


#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
    ${Teem_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../common
 )

include_directories( ${MODULE_INCLUDE_DIRECTORIES} )

set(MODULE_SRCS
    tractography.cc
    unscented_kalman_filter.cc
    filter_model.cc
    seed.cc
    ukffiber.cc
    NrrdData.cc
    vtk_writer.cc
    dwi_normalize.cc
    utilities.cc
    thread.cc
    QuadProg++_vnl.cc
    )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    teem
    ${CMAKE_THREAD_LIBS_INIT}
  )

#-----------------------------------------------------------------------------
enable_testing()
add_subdirectory(Testing)

#-----------------------------------------------------------------------------
if(UKFTractography_SuperBuild)
    SEMMacroBuildCLI(
        NAME ${MODULE_NAME}
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
        EXECUTABLE_ONLY
        )
endif()
