get_filename_component(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Input REALPATH)
get_filename_component(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../Data/Baseline REALPATH)
set(TESTING_RESULTS_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/Temporary")
file(MAKE_DIRECTORY ${TESTING_RESULTS_DIRECTORY}) ## Need to generate the directory for writing the temporary test outputs to

set(CLP UKFTractography)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

add_executable( ${CLP}Test CompareFibers.cc )
target_link_libraries( ${CLP}Test ${VTK_LIBRARIES} ${ITK_LIBRARIES})
add_dependencies(${CLP}Test ${CLP})
set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})


#
# remove output files before running tests
macro(RMTestFile CleanupTestName FileToRemove)
  add_test(NAME ${CleanupTestName}
    COMMAND ${CMAKE_COMMAND} -E remove ${FileToRemove}
    )
endmacro()




##############################################################################
# UKFTractography tests
##############################################################################


##############################################################################
# 1T no free water
# ----------------
set(testname ${CLP}_1T_Test)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/1T_fiber.vtk)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
  --dwiFile ${INPUT}/single_tensor.nhdr
  --maskFile ${INPUT}/mask.nhdr
  --tracts ${TESTING_RESULTS_DIRECTORY}/1T_fiber.vtk
  --seedsFile ${INPUT}/seed.nhdr
  --seedsPerVoxel 1
  --numTensor 1
  --numThreads 1
  --minBranchingAngle 0.0
  --maxBranchingAngle 0.0
  --recordNMSE
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)


set(testname ${CLP}_1T_Test_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/1T_fiber.vtk
  ${BASELINE}/1T_fiber.vtk
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_1T_Test)


##############################################################################
# 1T free water
# -------------
set(testname ${CLP}_1T_fw_Test)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/1T_fw_fiber.vtk)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
  --dwiFile ${INPUT}/single_tensor_fw.nhdr
  --maskFile ${INPUT}/mask.nhdr
  --tracts ${TESTING_RESULTS_DIRECTORY}/1T_fw_fiber.vtk
  --seedsFile ${INPUT}/seed.nhdr
  --seedsPerVoxel 1
  --numTensor 1
  --numThreads 1
  --minBranchingAngle 0.0
  --maxBranchingAngle 0.0
  --recordNMSE
  --freeWater
  --recordFreeWater
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)

set(testname ${CLP}_1T_fw_Test_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/1T_fw_fiber.vtk
  ${BASELINE}/1T_fw_fiber.vtk
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_1T_fw_Test)


##############################################################################
# 2T no freewater
# ---------------
set(testname ${CLP}_2T_Test)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/2T_fiber.vtk)

add_test(NAME ${testname}
    COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
    --dwiFile ${INPUT}/two_tensor.nhdr
    --maskFile ${INPUT}/mask.nhdr
    --tracts ${TESTING_RESULTS_DIRECTORY}/2T_fiber.vtk
 --seedsFile ${INPUT}/seed.nhdr
 --seedsPerVoxel 1
 --numTensor 2
 --numThreads 1
 --minBranchingAngle 0.0
 --maxBranchingAngle 0.0
 --recordNMSE
 --stoppingFA 0.10
 --stoppingThreshold 0.05
 --Qm 0.001
 --Ql 10
 --Rs 0.015
  )
set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)

set(testname ${CLP}_2T_Test_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fiber.vtk
  ${BASELINE}/2T_fiber.vtk
  )

set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_2T_Test)


##############################################################################
# 2T free water
# -------------
set(testname ${CLP}_2T_fw_Test)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
  -a
  --dwiFile ${INPUT}/two_tensor_fw.nhdr
  --maskFile ${INPUT}/mask.nhdr
  --tracts ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  --seedsFile ${INPUT}/seed.nhdr
  --seedsPerVoxel 1
  --numTensor 2
  --numThreads 1
  --minBranchingAngle 0.0
  --maxBranchingAngle 0.0
  --recordNMSE
  --freeWater
  --recordFreeWater
  --stoppingFA 0.1
  --stoppingThreshold 0.05
  --Qm 0.01
  --Ql 10
  --Rs 0.015
  )
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)

set(testname ${CLP}_2T_fw_Test_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  ${BASELINE}/2T_fw_fiber.vtk
  )

set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_2T_fw_Test)

set(testname ${CLP}_2T_fw_TestBinary)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtk)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
  --dwiFile ${INPUT}/two_tensor_fw.nhdr
  --maskFile ${INPUT}/mask.nhdr
  --tracts ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtk
  --seedsFile ${INPUT}/seed.nhdr
  --seedsPerVoxel 1
  --numTensor 2
  --numThreads 1
  --minBranchingAngle 0.0
  --maxBranchingAngle 0.0
  --recordNMSE
  --freeWater
  --recordFreeWater
  --stoppingFA 0.1
  --stoppingThreshold 0.05
  --Qm 0.01
  --Ql 10
  --Rs 0.015
  )
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)

set(testname ${CLP}_2T_fw_TestBinary_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtk
  ${BASELINE}/2T_fw_fiber.vtk
  )

set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_2T_fw_Test)
set_tests_properties(${testname} PROPERTIES DEPENDS ${CLP}_2T_fw_TestBinary)

set(testname ${CLP}_2T_fw_TestBinaryVTP)
RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtp)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
  --dwiFile ${INPUT}/two_tensor_fw.nhdr
  --maskFile ${INPUT}/mask.nhdr
  --tracts ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtp
  --seedsFile ${INPUT}/seed.nhdr
  --seedsPerVoxel 1
  --numTensor 2
  --numThreads 1
  --minBranchingAngle 0.0
  --maxBranchingAngle 0.0
  --recordNMSE
  --freeWater
  --recordFreeWater
  --stoppingFA 0.1
  --stoppingThreshold 0.05
  --Qm 0.01
  --Ql 10
  --Rs 0.015
  )
set_tests_properties(${testname} PROPERTIES DEPENDS ${testname}-cleanup)

set(testname ${CLP}_2T_fw_TestBinaryVTP_Compare)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber-binary.vtp
  ${BASELINE}/2T_fw_fiber.vtk
  )

set_property(TEST ${testname} PROPERTY LABELS ${CLP})
set_tests_properties(${testname} PROPERTIES DEPENDS "${CLP}_2T_fw_TestBinary;${CLP}_2T_fw_TestBinaryVTP")



##############################################################################
# Dummy Test as checkpoint to prevent races.
  add_test(NAME DUMMY_TEST
           COMMAND cmake -E sleep 1)
  set_tests_properties(DUMMY_TEST PROPERTIES DEPENDS
    "UKFTractography_1T_Test_Compare;UKFTractography_2T_fw_Test_Compare;UKFTractography_2T_fw_TestBinary_Compare;UKFTractography_2T_fw_TestBinaryVTP_Compare")
  set_tests_properties(DUMMY_TEST PROPERTIES RUN_SERIAL false)

##############################################################################


##############################################################################
# ConvertVTK tests
##############################################################################


####
# convert from legacy vtk ASCII to legacy vtk binary
####

set(testname ConvertVTK_LegacyText_ToLegacyBinary)
add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:ConvertVTK>
  -i ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  -o ${TESTING_RESULTS_DIRECTORY}/${testname}.vtk
  )

# depends on last of chain of 2T_fw_fiber tests; will depend on all of them
set_tests_properties(${testname} PROPERTIES DEPENDS "DUMMY_TEST;${testname}-cleanup")

add_test(NAME ${testname}Compare
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  ${TESTING_RESULTS_DIRECTORY}/${testname}.vtk
  )
set_property(TEST ${testname}Compare PROPERTY LABELS ConvertVTK)
set_tests_properties(${testname}Compare PROPERTIES DEPENDS ${testname})


####
# convert from legacy vtk ASCII to XML ASCII
####


set(testname ConvertVTK_LegacyText_ToXML)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:ConvertVTK>
  -a -u
  -i ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  -o ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp
  )

# depends on last of chain of 2T_fw_fiber tests; will depend on all of them
set_tests_properties(${testname} PROPERTIES DEPENDS "DUMMY_TEST;${testname}-cleanup")

add_test(NAME ${testname}Compare
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp
  )

set_property(TEST ${testname}Compare PROPERTY LABELS ConvertVTK)
set_tests_properties(${testname}Compare PROPERTIES DEPENDS ${testname})


####
# convert from legacy vtk ASCII to XML compressed binary
####

set(testname ConvertVTK_LegacyText_ToXMLBinaryCompressed)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:ConvertVTK>
  -i ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  -o ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp
  )

set_tests_properties(${testname} PROPERTIES DEPENDS "DUMMY_TEST;${testname}-cleanup")

add_test(NAME ${testname}Compare
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/2T_fw_fiber.vtk
  ${TESTING_RESULTS_DIRECTORY}/${testname}.vtp
  )

set_property(TEST ${testname}Compare PROPERTY LABELS ConvertVTK)
set_tests_properties(${testname}Compare PROPERTIES DEPENDS ${testname})

set(XMLBinaryCompressed ${testname})


####
# convert from XML compressed binary to legacy vtk ASCII
####

set(testname ConvertVTK_XMLBinaryCompressed_To_LegacyText)

RMTestFile(${testname}-cleanup ${TESTING_RESULTS_DIRECTORY}/${testname}.vtk)

add_test(NAME ${testname}
  COMMAND ${SEM_LAUNCH_COMMAND} $<TARGET_FILE:ConvertVTK>
  -a
  -i ${TESTING_RESULTS_DIRECTORY}/${XMLBinaryCompressed}.vtp
  -o ${TESTING_RESULTS_DIRECTORY}/${testname}.vtk
  )
# depends on last of chain of 2T_fw_fiber tests; will depend on all of them
set_tests_properties(${testname} PROPERTIES DEPENDS "DUMMY_TEST;${testname}-cleanup;${XMLBinaryCompressed}")

add_test(NAME ${testname}Compare
  COMMAND ${SEM_LAUNCH_COMMAND} ${CLP}Test
  ${TESTING_RESULTS_DIRECTORY}/${XMLBinaryCompressed}.vtp
  ${TESTING_RESULTS_DIRECTORY}/${testname}.vtk
  )

set_property(TEST ${testname}Compare PROPERTY LABELS ConvertVTK)
set_tests_properties(${testname}Compare PROPERTIES DEPENDS ${testname})

###############################################################################
# UKF Slicer CLI Test
###############################################################################
if (${UKFTractography_BUILD_SLICER_EXTENSION})
    if (MSVC) # need release subdir on Windows :|
      set(REL_VAR "/Release")
    else()
      set(REL_VAR "")
    endif()

  set(UKF_UNITTEST_LIB_PATHS "--additional-module-paths;${CMAKE_BINARY_DIR}/${Slicer_QTSCRIPTEDMODULES_LIB_DIR}${REL_VAR};${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}${REL_VAR};${CMAKE_BINARY_DIR}/${Slicer_CLIMODULES_LIB_DIR}${REL_VAR}")

    set(UKF_TEST_CLI_PYTHON "${CMAKE_CURRENT_SOURCE_DIR}/slicer_cli_test.py")

    set(testname UKF_Slicer_CLI)
    add_test(
      NAME ${testname}
      COMMAND ${Slicer_LAUNCHER_EXECUTABLE}
      ${UKF_UNITTEST_LIB_PATHS}
      --python-code "data_dir = '${INPUT}'; exec(open('${UKF_TEST_CLI_PYTHON}').read(), globals())"
      )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 180)

endif()
