# Works on Unix-like systems if the build folder within the root of source code and named 'build'

start=`date +%s`

SRC="../UKFTractography"

#logFile Name
logFile="log.txt"

# BINARY
BINARY='../build/UKFTractography-build/UKFTractography/bin/UKFTractography'

# VOLUME
dwi_path="$SRC/Data/Input/dwi.nhdr"

# MASK
mask_path="$SRC/Data/Input/dwi-mask.nhdr"

# SEEDS
seeds_path="$SRC/Data/Input/seeds_tc.nhdr"

# OUTPUT FIBER
output_path='./output_tractography.vtk'

eval $BINARY \
 --dwiFile $dwi_path \
 --maskFile $mask_path \
 --tracts $output_path \
 --seedsFile $seeds_path \
 --seedsPerVoxel 5 \
 --numTensor 2 | tee -a $logFile
 end=`date +%s`

runtime=$(python -c "print(${end} - ${start})")
echo "Output file name $output_path" | tee -a $logFile
echo "CPU Runtime was $runtime seconds"  | tee -a $logFile
