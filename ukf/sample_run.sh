#!/bin/bash

SRC="../ukf_tractography"

# BINARY
BINARY='bin/UKFTractography'

# VOLUME
diff_im_path="$SRC/data/dwi.nhdr"

# MASK
mask_path="$SRC/data/dwi-mask.nhdr"

# SEEDS
seeds_path="$SRC/data/seeds_tc.nhdr"

# OUTPUT FIBER
output_path='./fiber.vtk'

eval $BINARY \
 --dwiFile $diff_im_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --seedsFile $seeds_path \
 --minBranchingAngle 0.0 \
 --maxBranchingAngle 0.0 \
 --seedsPerVoxel 5 \
 --numTensor 2  \
 --simpleTensorModel \

diff $output_path $SRC/data/fiber_gold.vtk
