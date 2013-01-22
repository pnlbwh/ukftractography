#-------------------- UKF Tractography Testing ---------------------------------
# Running four different synthetic data cases (1 tensor, 1 tensor and free water,
# 2 tensor, 2 tensor and free water) and comparing the output to precalculated
# fibers. The last test is deleting the output fiber files again in case the
# test has to be rerun.
#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Setting up the Binary for fiber comparison 
#-----------------------------------------------------------------------------

## Build Macros
include(${PROJECT_SOURCE_DIR}/BuildTools/BuildMacros.cmake)


## This variable is only available in CMAKE 2.8.3 
SETIFEMPTY(CMAKE_CURRENT_LIST_DIR ${PROJECT_SOURCE_DIR}/testing)

## The CLP Name
set(CLP UKFTractography)

## Set the directories for testing
set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/../data)
set(OUTPUT_DIR ${PROJECT_SOURCE_DIR}/testing)
set(TEST_CODE_DIR ${CMAKE_CURRENT_LIST_DIR})

## Build the code for testing
set(Compare_SOURCE 
    ${TEST_CODE_DIR}/CompareFibers.cc 
    ${TEST_CODE_DIR}/vtkReader.cc
    )

#include_directories(${PROJECT_SOURCE_DIR}/ukf)
add_executable(${CLP}Compare ${Compare_SOURCE})
add_dependencies(${CLP}Compare ${CLP})

#if (Slicer_VERSION EQUAL 4)
set(Slicer_EXE ${Slicer_HOME}/Slicer)
#endif()

#if (Slicer_VERSION EQUAL 3)
  #set(Slicer_EXE ${Slicer3_HOME}/Slicer3)
#endif()

## Add the correct executables
if(NOT Slicer_VERSION) 
    set(RUNEXEC 
        ${PROJECT_BINARY_DIR}/bin/${CLP}
        )
else()
    set(RUNEXEC 
        ${Slicer_EXE} --launch ${CLP}
        ) 
endif()

set(COMPEXEC ${PROJECT_BINARY_DIR}/testing/${CLP}Compare)

#-----------------------------------------------------------------------------
# Run the actual tests
#-----------------------------------------------------------------------------

######## 1T no free water #########################
add_test(${PROJECT_NAME}_1T_Test
 ${RUNEXEC}
 --dwiFile ${TEST_DATA_DIR}/single_tensor.nhdr
 --maskFile ${TEST_DATA_DIR}/mask.nhdr
 --tracts ${OUTPUT_DIR}/1T_fiber.vtk
 --seedsFile ${TEST_DATA_DIR}/seed.nhdr
 --seedsPerVoxel 1
 --numTensor 1
 --numThreads 1
 --minBranchingAngle 0.0
 --maxBranchingAngle 0.0
 --simpleTensorModel
 --recordNMSE
)
add_test(${PROJECT_NAME}_1T_Test_Compare
 ${COMPEXEC}
 ${OUTPUT_DIR}/1T_fiber.vtk
 ${TEST_DATA_DIR}/1T_fiber.vtk
)

######## 1T free water #############################
add_test(${PROJECT_NAME}_1T_fw_Test
 ${RUNEXEC}
 --dwiFile ${TEST_DATA_DIR}/single_tensor_fw.nhdr
 --maskFile ${TEST_DATA_DIR}/mask.nhdr
 --tracts ${OUTPUT_DIR}/1T_fw_fiber.vtk
 --seedsFile ${TEST_DATA_DIR}/seed.nhdr
 --seedsPerVoxel 1
 --numTensor 1
 --numThreads 1
 --minBranchingAngle 0.0
 --maxBranchingAngle 0.0
 --simpleTensorModel
 --recordNMSE
 --freeWater
 --recordFreeWater
)

add_test(${PROJECT_NAME}_1T_fw_Test_Compare
 ${COMPEXEC}
 ${OUTPUT_DIR}/1T_fw_fiber.vtk
 ${TEST_DATA_DIR}/1T_fw_fiber.vtk
)

######## 2T no freewater #############################
add_test(${PROJECT_NAME}_2T_Test
 ${RUNEXEC}
 --dwiFile ${TEST_DATA_DIR}/two_tensor.nhdr
 --maskFile ${TEST_DATA_DIR}/mask.nhdr
 --tracts ${OUTPUT_DIR}/2T_fiber.vtk
 --seedsFile ${TEST_DATA_DIR}/seed.nhdr
 --seedsPerVoxel 1
 --numTensor 2
 --numThreads 1
 --minBranchingAngle 0.0
 --maxBranchingAngle 0.0
 --simpleTensorModel
 --recordNMSE
 --minFA 0.10
 --minGA 0.05
 --Qm 0.001
 --Ql 10
 --Rs 0.015
)

add_test(${PROJECT_NAME}_2T_Test_Compare
 ${COMPEXEC}
 ${OUTPUT_DIR}/2T_fiber.vtk
 ${TEST_DATA_DIR}/2T_fiber.vtk
)

####### 2T free water #################################
add_test(${PROJECT_NAME}_2T_fw_Test
 ${RUNEXEC}
 --dwiFile ${TEST_DATA_DIR}/two_tensor_fw.nhdr
 --maskFile ${TEST_DATA_DIR}/mask.nhdr
 --tracts ${OUTPUT_DIR}/2T_fw_fiber.vtk
 --seedsFile ${TEST_DATA_DIR}/seed.nhdr
 --seedsPerVoxel 1
 --numTensor 2
 --numThreads 1
 --minBranchingAngle 0.0
 --maxBranchingAngle 0.0
 --simpleTensorModel
 --recordNMSE
 --freeWater
 --recordFreeWater
 --minFA 0.1
 --minGA 0.05
 --Qm 0.01
 --Ql 10
 --Rs 0.015
)

add_test(${PROJECT_NAME}_2T_fw_Test_Compare
 ${COMPEXEC}
 ${OUTPUT_DIR}/2T_fw_fiber.vtk
 ${TEST_DATA_DIR}/2T_fw_fiber.vtk
)

# Clean up the generated files
add_test(${PROJECT_NAME}_cleanup
 rm
 ${OUTPUT_DIR}/1T_fiber.vtk
 ${OUTPUT_DIR}/1T_fw_fiber.vtk
 ${OUTPUT_DIR}/2T_fiber.vtk
 ${OUTPUT_DIR}/2T_fw_fiber.vtk
)


