#-----------------------------------------------------------------------------
set(MODULE_NAME UKFTractography)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  ${Teem_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  )

set(MODULE_SRCS
  BRAINSThreadControl.cxx
  tractography.cc
  unscented_kalman_filter.cc
  seed.cc
  ukffiber.cc
  NrrdData.cc
  vtk_writer.cc
  dwi_normalize.cc
  utilities.cc
  thread.cc
  QuadProg++_Eigen.cc
  filter_model.cc
  filter_Full1T.cc
  filter_Full1T_FW.cc
  filter_Full2T.cc
  filter_Full2T_FW.cc
  filter_Full3T.cc
  filter_NODDI1F.cc
  filter_NODDI2F.cc
  filter_Simple1T.cc
  filter_Simple1T_FW.cc
  filter_Simple2T.cc
  filter_Simple2T_FW.cc
  filter_Simple3T.cc
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${TEEM_LIB} ${ZLIB_LIBRARIES}
  ${VTK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
add_subdirectory(Testing)

message( STATUS "XXXX ${SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION}")
#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ConvertVTK
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  EXECUTABLE_ONLY
  )
