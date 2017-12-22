/**
 * \file ISignalData.h
 * \brief Contains interface for signal data. Currently implemented by NrrdData.* using teem
 * \todo Reimplement NrrdData using ITK or VTK
 * \todo Rethink what really belongs into the interface, and what is implementation specific
*/

#ifndef ISIGNALDATA_H_
#define ISIGNALDATA_H_

#include <iostream>
#include <string>
#include <vector>
#include <teem/nrrd.h>
#include "linalg.h"

/**
 * \class ISignalData
 * \brief Generic Interface for signal data
 *
 * Contains the function definitions for working with the diffusion and label data
 * on which the rest of the code relies. It's current implementation is NrrdData.*
*/
class ISignalData
{
public:

  /**
   * Constructor
   * \param[in] sigma_signal the interpolation 'factor' for the signal
   * \param[in] sigma_mask the interpolation 'factor' for the mask
  */
  ISignalData(ukfPrecisionType sigma_signal, ukfPrecisionType sigma_mask)
    : _sigma_signal(sigma_signal), _sigma_mask(sigma_mask)
  {

  }

  /** Deconstructor */
  virtual ~ISignalData()
  {
  }

  /** Gets the signal values at a specified position. */
  virtual void Interp3Signal(const vec3_t& pos, ukfVectorType & signal) const = 0;

  /** Checks if a certian position is still within the brain mask. */
  virtual ukfPrecisionType Interp3ScalarMask(const vec3_t& pos) const = 0;

  /** Checks if a certian position is still within the brain mask. */
  virtual ukfPrecisionType ScalarMaskValue(const vec3_t& pos) const = 0;

  /** Get all the seed points. */
  virtual void GetSeeds(const std::vector<int>& labels, stdVec_t& seeds) const = 0;

  /** Returns the gradients. */
  virtual const stdVec_t & gradients() const = 0;

  /** Returns the vector of b values */
  virtual const ukfVectorType & GetBValues() const = 0;

  /** Return the size of the signal vector */
  virtual int GetSignalDimension() const = 0;

  /**
    * \brief Load all Data
    * \param[in] data_file The path of the diffusion image
    * \param[in] seed_file The path of the seeds, a binary label map containing the starting points
    * \param[in] mask_file The path of the non-optional brain mask
    * \param[in] normalizedDWIData If set to 'true', the data will not be normalized
    * \param[in] outputNormalizedDWIData If set to 'true' the result of the normalization will be saved
    *
    * Loads all the data necessary to perform tractography
  */
  virtual bool LoadData(const std::string& data_file, const std::string& seed_file, const std::string& mask_file,
                        const bool normalizedDWIData, const bool outputNormalizedDWIData) = 0;

  /** Returns the dimensions of the image */
  virtual vec3_t dim() const = 0;

  /** Returns the voxel spacing */
  vec3_t voxel() const
  {
    return _voxel;
  }

  /** Returns the ijk-to-RAS matrix */
  const ukfMatrixType i2r() const
  {
    return _i2r;
  }

  /** Returns the RAS-to-ijk matrix */
  const ukfMatrixType r2i() const
  {
    return _r2i;
  }

protected:

  /** sigma for gaussian interpolation of signal */
  const ukfPrecisionType _sigma_signal;
  /** sigma for gaussian interpolation of mask */
  const ukfPrecisionType _sigma_mask;

  /** voxel size */
  vec3_t _voxel;

  /** matrix for RAS to ijk conversion */
  ukfMatrixType _r2i;

  /** matrix for ijk to RAS conversion */
  ukfMatrixType _i2r;
};

#endif // ISIGNALDATA_H
