#!/bin/bash

## NEEDS TO BE RUN FROM THE test_code FOLDER!

# BINARY
#ukf_path="../../ukf-build/bin/UKFTractography"

# MASK
mask_path='Data/Input/mask.nhdr'

# SEEDS
seeds_path='Data/Input/seed.nhdr'

## 1T
diff_im_path='Data/Input/single_tensor.nhdr'
output_path='Data/Baseline/1T_fiber.vtk'
my_command="$ukf_path --dwiFile $diff_im_path --maskFile $mask_path --tracts $output_path --seedsFile $seeds_path
--seedsPerVoxel 1  --numTensor 1  --numThreads 1 --minBranchingAngle 0.0 --maxBranchingAngle 0.0 --recordNMSE"
eval $my_command

## 1T - fw
diff_im_path='Data/Input/single_tensor_fw.nhdr'
output_path='Data/Baseline/1T_fw_fiber.vtk'
my_command="$ukf_path --dwiFile $diff_im_path --maskFile $mask_path --tracts $output_path --seedsFile $seeds_path
--seedsPerVoxel 1  --numTensor 1  --numThreads 1 --minBranchingAngle 0.0 --maxBranchingAngle 0.0 --freeWater --recordFreeWater --recordNMSE"
eval $my_command

## 2T
diff_im_path='Data/Input/two_tensor.nhdr'
output_path='Data/Baseline/2T_fiber.vtk'
my_command="$ukf_path --dwiFile $diff_im_path --maskFile $mask_path --tracts $output_path --seedsFile $seeds_path
--seedsPerVoxel 1  --numTensor 2  --numThreads 1 --minBranchingAngle 0.0 --maxBranchingAngle 0.0 --stoppingFA 0.10 --stoppingThreshold 0.05 --Qm 0.001 --Ql 10 --Rs 0.015 --recordNMSE"
eval $my_command
#
## 2T - fw
diff_im_path='Data/Input/two_tensor_fw.nhdr'
output_path='Data/Baseline/2T_fw_fiber.vtk'
my_command="$ukf_path --dwiFile $diff_im_path --maskFile $mask_path --tracts $output_path --seedsFile $seeds_path
--seedsPerVoxel 1  --numTensor 2  --numThreads 1 --minBranchingAngle 0.0 --maxBranchingAngle 0.0 --freeWater --recordFreeWater --recordNMSE --stoppingFA 0.10 --stoppingThreshold 0.05 --Qm 0.01 --Ql 10 --Rs 0.015"
eval $my_command

