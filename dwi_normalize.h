/**
 * \file dwi_normalize.h
 * \brief Contains function for preprocessing of DW Data
 * \author Yinpeng Li (mousquetaires@unc.edu)
*/

#ifndef DWI_NORMALIZE_H_
#define DWI_NORMALIZE_H_

#include <teem/nrrd.h>

/**
 * \brief Preprocess diffusion data for tractography
 *
 * This function performs the necessary preprocessings for the ukf tractography
 * It computes baseline image by avaraging, and divides out baseline image from signal
 * It also adjusts the axis order in NRRD, and transforms the world coordinate system to RAS if necessary
 *
 * This function complies with the NRRD convention, namely the parameters have to be initialized properly
 * before invocation
*/

void dwiNormalize(const Nrrd *raw, Nrrd *&normalized) ;

const int DATA_DIMENSION = 4 ;

#endif
