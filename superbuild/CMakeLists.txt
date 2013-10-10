# ---------- Superbuild system for UKFTractography -------------------------------------------------------
# SuperBuild script for the Unscented Kalman Filter Project. It takes care of all dependancies and
# then invokes the main project
# --------------------------------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.4)
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

project(UKFTractographySuperBuild)

include(ExternalProject)

option(BUILD_SHARED_LIBS "Build with shared libraries." OFF)

#-----------------------------------------------------------------------------
# Settings and Options
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

string(TOUPPER CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} BUILD_C_FLAGS)
string(TOUPPER CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} BUILD_CXX_FLAGS)


set(COMMON_CMAKE_ARGS
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  )

#-----------------------------------------------------------------------------
# Get External Dependacies
#-----------------------------------------------------------------------------
ExternalProject_add(ITK
  GIT_REPOSITORY git://itk.org/ITK.git
  GIT_TAG v4.0b01
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DITK_LEGACY_REMOVE:BOOL=ON
  -DITK_BUILD_ALL_MODULES:BOOL=ON
  -DINSTALL_WRAP_ITK_COMPATIBILITY:BOOL=OFF
  )
set(ITK_DIR ${CMAKE_CURRENT_BINARY_DIR}/ITK-prefix/src/ITK-build)

ExternalProject_add(SlicerExecutionModel
  GIT_REPOSITORY https://github.com/Slicer/SlicerExecutionModel.git
  GIT_TAG "origin/master"
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS ${COMMON_CMAKE_ARGS}
  -DITK_DIR:PATH=${ITK_DIR}
  DEPENDS ITK
  )
set(SlicerExecutionModel_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/SlicerExecutionModel-prefix/src/SlicerExecutionModel-build)
set(GenerateCLP_DIR ${SlicerExecutionModel_DIR}/GenerateCLP)

# the 'find_package' files are not present in the build directory so you have
# to install teem somewhere for find_package to work
set(TEEM_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/teem-prefix/src/teem-install)
ExternalProject_add(teem
  SVN_REPOSITORY "http://svn.code.sf.net/p/teem/code/teem/tags/1.10.0"
  UPDATE_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DCMAKE_INSTALL_PREFIX:PATH=${TEEM_INSTALL_PREFIX}
  -DTeem_USE_LIB_INSTALL_SUBDIR:BOOL=OFF
  -DTeem_PTHREAD:BOOL=OFF
  -DTeem_BZIP2:BOOL=OFF
  -DTeem_ZLIB:BOOL=ON
  -DTeem_PNG:BOOL=OFF
  # -DTeem_VTK_MANGLE:BOOL=ON
  )
set(Teem_DIR ${TEEM_INSTALL_PREFIX}/lib/Teem-1.10.0)
#set(Teem_DIR ${TEEM_INSTALL_PREFIX}/lib/)
#set(teem_INCLUDE_DIR ${TEEM_INSTALL_PREFIX}/include)

#set(BOOST_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/boost-prefix/src/boost-install)
#ExternalProject_add(boost
  #URL http://downloads.sourceforge.net/project/boost/boost/1.41.0/boost_1_41_0.tar.gz
  #URL_MD5 ddb04c9b9164a36a91dcf36819666dd9
  #BUILD_IN_SOURCE 1
  #CONFIGURE_COMMAND
  #<SOURCE_DIR>/bootstrap.sh
  #--prefix=${BOOST_INSTALL_PREFIX}
  #link=static
  #runtime-link=static
  #BUILD_COMMAND "" # <SOURCE_DIR>/bjam link=static runtime-link=static
  #INSTALL_COMMAND <SOURCE_DIR>/bjam link=static runtime-link=static install
#)
#set(Boost_INCLUDE_DIRS ${BOOST_INSTALL_PREFIX}/include)
#set(Boost_LIBRARY_DIRS ${BOOST_INSTALL_PREFIX}/lib)

#-----------------------------------------------------------------------------
# Invoke the main project
#-----------------------------------------------------------------------------
ExternalProject_add(UKFCommonLib
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../common/
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DUKFTractography_SuperBuild::bool=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DTeem_DIR:PATH=${Teem_DIR}
  -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ITK teem SlicerExecutionModel
)

ExternalProject_add(UKFTractography
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../ukf/
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DUKFTractography_SuperBuild::bool=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DTeem_DIR:PATH=${Teem_DIR}
  -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ITK teem SlicerExecutionModel UKFCommonLib
)

ExternalProject_add(vtk2mask
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../vtk2mask/
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DUKFTractography_SuperBuild::bool=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ITK SlicerExecutionModel UKFCommonLib
)

ExternalProject_add(vtkFilter
  SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../vtkFilter/
  CMAKE_ARGS
  ${COMMON_CMAKE_ARGS}
  -DUKFTractography_SuperBuild::bool=ON
  -DITK_DIR:PATH=${ITK_DIR}
  -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
  -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ITK SlicerExecutionModel UKFCommonLib
)
