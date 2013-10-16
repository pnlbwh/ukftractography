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
#include <vnl/vnl_matrix.h>

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
  ISignalData(double sigma_signal, double sigma_mask)
    : _sigma_signal(sigma_signal), _sigma_mask(sigma_mask)
  {

  }

  /** Deconstructor */
  virtual ~ISignalData()
  {
  }

  /** Gets the signal values at a specified position. */
  virtual void Interp3Signal(const vec_t& pos, std::vector<double>& signal) const = 0;

  /** Checks if a certian position is still within the brain mask. */
  virtual double Interp3ScalarMask(const vec_t& pos) const = 0;

  /** Get all the seed points. */
  virtual void GetSeeds(const std::vector<int>& labels, std::vector<vec_t>& seeds) const = 0;

  /** Returns the gradients. */
  virtual const std::vector<vec_t> & gradients() const = 0;

  /** Returns the vector of b values */
  virtual const std::vector<double> & GetBValues() const = 0;

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

  /**
   * Loading the Signal
   * \todo this shouldn't be in the interface, because it is only used in the implementation, make it a private function of NrrdData
  */
  virtual bool LoadSignal(const std::string& data_file, const bool normalizedDWIData) = 0;

  /** Returns the voxel spacing */
  vec_t voxel() const
  {
    return _voxel;
  }

  /** Returns the dimensions of the image */
  virtual vec_t dim() const = 0;

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
  const double _sigma_signal;
  /** sigma for gaussian interpolation of mask */
  const double _sigma_mask;

  /** voxel size */
  vec_t _voxel;

  /** matrix for RAS to ijk conversion */
  ukfMatrixType _r2i;

  /** matrix for ijk to RAS conversion */
  ukfMatrixType _i2r;
};

#endif // ISIGNALDATA_H
