#-----------------------------------------------------------------------------
PROJECT(UKFBase)

#-----------------------------------------------------------------------------
# Generate cpp with git commit hash
set(UKF_GITVER_FILE "${CMAKE_CURRENT_BINARY_DIR}/git_version.cc")
set(UKF_GITVER_TEMP "${CMAKE_CURRENT_BINARY_DIR}/git_version.tmp")

# always run this target to check for changes.
add_custom_target(_GEN_GITVER ALL
    # write git hash to temp file
    COMMAND ${CMAKE_COMMAND} -E echo         "#include \"git_version.h\""                            > "${UKF_GITVER_TEMP}"
    COMMAND ${CMAKE_COMMAND} -E echo_append  "#define _UKF_GIT_HASH "                               >> "${UKF_GITVER_TEMP}"
    COMMAND                                  ${GIT_EXECUTABLE} rev-parse HEAD                       >> "${UKF_GITVER_TEMP}"
    COMMAND ${CMAKE_COMMAND} -E echo_append  "const char UKF_GIT_HASH[] = _ukf_str(_UKF_GIT_HASH);" >> "${UKF_GITVER_TEMP}"

    # copy to the real file if different
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${UKF_GITVER_TEMP}" "${UKF_GITVER_FILE}"

    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating UKF_GIT_HASH"
    VERBATIM
    )

# dummy output so we can depend on it in the external macro, because CMake.
add_custom_command(OUTPUT ${UKF_GITVER_FILE}
                   DEPENDS _GEN_GITVER)

#-----------------------------------------------------------------------------
# SlicerExecutionModel
find_package(SlicerExecutionModel REQUIRED)
include(${GenerateCLP_USE_FILE})

set(UKFCLI_CXX "${UKF_SRC}/UKFTractography/UKFTractography.cxx")
# TODO remove
generateclp(UKFCLI_SOURCE "${UKF_SRC}/UKFTractography/UKFTractography.xml")

#-----------------------------------------------------------------------------
set(MODULE_SRCS
  ../common/BRAINSThreadControl.cxx
  utilities.cc
  cli.cc
  tractography.cc
  unscented_kalman_filter.cc
  seed.cc
  ukffiber.cc
  NrrdData.cc
  vtk_writer.cc
  dwi_normalize.cc
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
  "${UKF_GITVER_FILE}"
  UKFTractographyCLP.h
  )

set(MODULE_INCLUDE_DIRECTORIES
  ${Teem_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  ${CMAKE_CURRENT_SOURCE_DIR}/../UKFTractography
  ${CMAKE_CURRENT_BINARY_DIR}
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${TEEM_LIB} ${ZLIB_LIBRARIES}
  ${VTK_LIBRARIES}
  ${SlicerExecutionModel_LIBRARIES}
  )

add_library(${PROJECT_NAME}
            SHARED
            ${MODULE_SRCS})

# Add ordering dependency on Eigen external project target.
add_dependencies(UKFBase EP_Eigen)

target_include_directories(${PROJECT_NAME} PUBLIC
                           ${MODULE_INCLUDE_DIRECTORIES})

target_link_libraries(${PROJECT_NAME} ${MODULE_TARGET_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
    )
