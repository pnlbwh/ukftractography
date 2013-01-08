#!/bin/bash

# BINARY
BINARY='bin/UKFTractography'

# VOLUME
diff_im_path='/projects/schiz/3Tdata/case01045/diff/01045-dwi-filt-Ed.nhdr'

# MASK
mask_path='/projects/schiz/3Tdata/case01045/diff/Tensor_mask-01045-dwi-filt-Ed_AvGradient-edited.nhdr'

# SEEDS
seeds_path='/home/malcolm/src/seeds_tc.nhdr'

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

