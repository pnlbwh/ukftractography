#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)


#-----------------------------------------------------------------------------
# Where should the superbuild source files be downloaded to?
# By keeping this outside of the build tree, you can share one
# set of external source trees for multiple build trees
#-----------------------------------------------------------------------------
set( SOURCE_DOWNLOAD_CACHE ${CMAKE_CURRENT_BINARY_DIR}/src CACHE PATH
    "The path for downloading external source directories" )
mark_as_advanced( SOURCE_DOWNLOAD_CACHE )

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------
include(ExternalProject)

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()

#-----------------------------------------------------------------------------
set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")
set_property(CACHE EXTERNAL_PROJECT_BUILD_TYPE PROPERTY
  STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_BOOST "Build using an externally defined version of BOOST" OFF)

#------------------------------------------------------------------------------
set(SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY bin)
set(SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY lib)
set(SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY lib)
set(SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION bin)
set(SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION lib)
set(SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION lib)

mark_as_superbuild(
  VARS
    SlicerExecutionModel_INSTALL_BIN_DIR:STRING
    SlicerExecutionModel_INSTALL_LIB_DIR:STRING
    SlicerExecutionModel_INSTALL_NO_DEVELOPMENT
  PROJECTS SlicerExecutionModel
  )

#------------------------------------------------------------------------------
# ${PRIMARY_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------

## for i in SuperBuild/*; do  echo $i |sed 's/.*External_\([a-zA-Z]*\).*/\1/g'|fgrep -v cmake|fgrep -v Template; done|sort -u
set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES
  Boost
  SlicerExecutionModel
  ITK
  Eigen
  VTK
  teem
  )

#-----------------------------------------------------------------------------
# Common external projects CMake variables
#-----------------------------------------------------------------------------
mark_as_superbuild(
  VARS
    MAKECOMMAND:STRING
    CMAKE_SKIP_RPATH:BOOL
    BUILD_SHARED_LIBS:BOOL
    CMAKE_CXX_COMPILER:PATH
    CMAKE_CXX_FLAGS:STRING
    CMAKE_CXX_FLAGS_DEBUG:STRING
    CMAKE_CXX_FLAGS_MINSIZEREL:STRING
    CMAKE_CXX_FLAGS_RELEASE:STRING
    CMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_C_COMPILER:PATH
    CMAKE_C_FLAGS:STRING
    CMAKE_C_FLAGS_DEBUG:STRING
    CMAKE_C_FLAGS_MINSIZEREL:STRING
    CMAKE_C_FLAGS_RELEASE:STRING
    CMAKE_C_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_EXE_LINKER_FLAGS:STRING
    CMAKE_EXE_LINKER_FLAGS_DEBUG:STRING
    CMAKE_EXE_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_EXE_LINKER_FLAGS_RELEASE:STRING
    CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_MODULE_LINKER_FLAGS:STRING
    CMAKE_MODULE_LINKER_FLAGS_DEBUG:STRING
    CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_MODULE_LINKER_FLAGS_RELEASE:STRING
    CMAKE_MODULE_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_SHARED_LINKER_FLAGS:STRING
    CMAKE_SHARED_LINKER_FLAGS_DEBUG:STRING
    CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL:STRING
    CMAKE_SHARED_LINKER_FLAGS_RELEASE:STRING
    CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO:STRING
    CMAKE_GENERATOR:STRING
    CMAKE_EXTRA_GENERATOR:STRING
    CMAKE_EXPORT_COMPILE_COMMANDS:BOOL
    CMAKE_INSTALL_PREFIX:PATH
    CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
    CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
    CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
    CMAKE_INSTALL_RUNTIME_DESTINATION:PATH
    CMAKE_INSTALL_LIBRARY_DESTINATION:PATH
    CMAKE_INSTALL_ARCHIVE_DESTINATION:PATH
    CTEST_NEW_FORMAT:BOOL
    MEMORYCHECK_COMMAND_OPTIONS:STRING
    MEMORYCHECK_COMMAND:PATH
    CMAKE_SHARED_LINKER_FLAGS:STRING
    CMAKE_EXE_LINKER_FLAGS:STRING
    CMAKE_MODULE_LINKER_FLAGS:STRING
    SITE:STRING
    BUILDNAME:STRING
    ${PROJECT_NAME}_BUILD_DICOM_SUPPORT:BOOL
    PYTHON_EXECUTABLE:FILEPATH
    PYTHON_INCLUDE_DIR:PATH
    PYTHON_LIBRARY:FILEPATH
    BOOST_ROOT:PATH
    BOOST_INCLUDE_DIR:PATH
    Boost_LIBRARY_DIR:PATH
    SlicerExecutionModel_DIR:PATH
    SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
    UKFTractography_USE_GIT_PROTOCOL:BOOL
  ALL_PROJECTS
  )

if(${PRIMARY_PROJECT_NAME}_USE_QT)
  mark_as_superbuild(
    VARS
      ${PRIMARY_PROJECT_NAME}_USE_QT:BOOL
      QT_QMAKE_EXECUTABLE:PATH
      QT_MOC_EXECUTABLE:PATH
      QT_UIC_EXECUTABLE:PATH
    ALL_PROJECTS
    )
endif()
mark_as_superbuild(${PRIMARY_PROJECT_NAME}_USE_QT)

#-----------------------------------------------------------------------------
# Set CMake OSX variable to pass down the external projects
#-----------------------------------------------------------------------------
if(APPLE)
  mark_as_superbuild(
    VARS
      CMAKE_OSX_ARCHITECTURES:STRING
      CMAKE_OSX_SYSROOT:PATH
      CMAKE_OSX_DEPLOYMENT_TARGET:STRING
    ALL_PROJECTS
    )
endif()

set(extProjName ${PRIMARY_PROJECT_NAME})
set(proj        ${PRIMARY_PROJECT_NAME})
ExternalProject_Include_Dependencies(${proj} DEPENDS_VAR ${PRIMARY_PROJECT_NAME}_DEPENDENCIES)

#-----------------------------------------------------------------------------
# Add external project CMake args
#-----------------------------------------------------------------------------

mark_as_superbuild(
  VARS
    BUILD_EXAMPLES:BOOL
    BUILD_TESTING:BOOL
  )

#-----------------------------------------------------------------------------
#
# By default we want to build ${PROJECT_NAME} stuff using the CMAKE_BUILD_TYPE of
# the top level build, but build the support libraries in Release.
# So make a list of option that will be passed only to all the prerequisite libraries.
#
set(COMMON_EXTERNAL_PROJECT_ARGS)
if(NOT CMAKE_CONFIGURATION_TYPES)
  list(APPEND COMMON_EXTERNAL_PROJECT_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${EXTERNAL_PROJECT_BUILD_TYPE}
    )
endif()

#-----------------------------------------------------------------------------
# CTestCustom
#-----------------------------------------------------------------------------
if(BUILD_TESTING AND NOT ${PRIMARY_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  configure_file(
    CMake/CTestCustom.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/CTestCustom.cmake
    @ONLY)
endif()

#------------------------------------------------------------------------------
# Configure and build ${PROJECT_NAME}
#------------------------------------------------------------------------------
set(proj ${PRIMARY_PROJECT_NAME})
ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DEPENDS ${${PRIMARY_PROJECT_NAME}_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${PRIMARY_PROJECT_NAME}-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    --no-warn-unused-cli    # HACK Only expected variables should be passed down.
    -D${PRIMARY_PROJECT_NAME}_SUPERBUILD:BOOL=OFF    #NOTE: VERY IMPORTANT reprocess top level CMakeList.txt
    -DRUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DLIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    -DARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
    -DINSTALL_RUNTIME_DESTINATION:PATH=${CMAKE_INSTALL_RUNTIME_DESTINATION}
    -DINSTALL_LIBRARY_DESTINATION:PATH=${CMAKE_INSTALL_LIBRARY_DESTINATION}
    -DINSTALL_ARCHIVE_DESTINATION:PATH=${CMAKE_INSTALL_ARCHIVE_DESTINATION}
  INSTALL_COMMAND ""
  )

# This custom external project step forces the build and later
# steps to run whenever a top level build is done...
ExternalProject_Add_Step(${proj} forcebuild
  COMMAND ${CMAKE_COMMAND} -E remove
    ${CMAKE_CURRENT_BINARY_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
  COMMENT "Forcing build step for '${proj}'"
  DEPENDEES build
  ALWAYS 1
  )

