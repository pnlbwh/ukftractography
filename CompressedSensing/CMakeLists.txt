
#-----------------------------------------------------------------------------
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
message("Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")

#-----------------------------------------------------------------------------
set(MODULE_NAME CompressedSensing)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  ${Teem_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(MODULE_SRCS
  CompressedSensing.cxx
  SphericalHarmonics.cxx
  DefaultGradients.cxx
  DoCSEstimate.cxx
  icosahedron3.cxx
  BuildRidges.cxx
  MultiSample.cxx
  BuildSensor.cxx
  SphHarm.cxx
  NrrdFile.cxx
  BPDN_homotopy.cxx
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${TEEM_LIB} ${ZLIB_LIBRARIES}
  ${VTK_LIBRARIES}
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
