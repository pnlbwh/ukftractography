
#-----------------------------------------------------------------------------
enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

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
# Git protocole option
#-----------------------------------------------------------------------------
option(${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
set(git_protocol "git")
if(NOT ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

find_package(Git REQUIRED)

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------
include(ExternalProject)

# With CMake 2.8.9 or later, the UPDATE_COMMAND is required for updates to occur.
# For earlier versions, we nullify the update state to prevent updates and
# undesirable rebuild.
option(FORCE_EXTERNAL_BUILDS "Force rebuilding of external project (if they are updated)" ON)
if(CMAKE_VERSION VERSION_LESS 2.8.9 OR NOT FORCE_EXTERNAL_BUILDS)
  set(cmakeversion_external_update UPDATE_COMMAND)
  set(cmakeversion_external_update_value "" )
else()
  set(cmakeversion_external_update LOG_UPDATE )
  set(cmakeversion_external_update_value 1)
endif()

#-----------------------------------------------------------------------------
set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)

#------------------------------------------------------------------------------
set(SlicerExecutionModel_INSTALL_BIN_DIR bin)
set(SlicerExecutionModel_INSTALL_LIB_DIR lib)
set(SlicerExecutionModel_INSTALL_NO_DEVELOPMENT 1)
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
    SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_LIBRARY_DESTINATION:PATH
    SlicerExecutionModel_DEFAULT_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  PROJECTS SlicerExecutionModel
  )

#------------------------------------------------------------------------------
# ${PRIMARY_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------
set(ITK_EXTERNAL_NAME ITKv${ITK_VERSION_MAJOR})

## for i in SuperBuild/*; do  echo $i |sed 's/.*External_\([a-zA-Z]*\).*/\1/g'|fgrep -v cmake|fgrep -v Template; done|sort -u
set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES
  SlicerExecutionModel
  ${ITK_EXTERNAL_NAME}
  VTK
  teem
  Eigen
  Boost
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
    CMAKE_INSTALL_PREFIX:PATH
    CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
    CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
    CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
    CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
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
  ALL_PROJECTS
  )

if(${PRIMARY_PROJECT_NAME}_USE_QT)
  mark_as_superbuild(VARS QT_QMAKE_EXECUTABLE QT_MOC_EXECUTABLE QT_UIC_EXECUTABLE)
endif()
mark_as_superbuild(${PRIMARY_PROJECT_NAME}_USE_QT)

#-----------------------------------------------------------------------------
# Set CMake OSX variable to pass down the external projects
#-----------------------------------------------------------------------------
if(APPLE)
  mark_as_superbuild(
    VARS CMAKE_OSX_ARCHITECTURES:STRING CMAKE_OSX_SYSROOT:PATH CMAKE_OSX_DEPLOYMENT_TARGET:STRING
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
    ITK_VERSION_MAJOR:STRING
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
  CMAKE_ARGS
    --no-warn-unused-cli    # HACK Only expected variables should be passed down.
    -D${PRIMARY_PROJECT_NAME}_SUPERBUILD:BOOL=OFF    #NOTE: VERY IMPORTANT reprocess top level CMakeList.txt
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

