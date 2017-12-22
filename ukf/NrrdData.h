/**
 * \file NrrdData.h
 * \brief Contains an implementation of ISignalData using the teem library.
*/

#ifndef NRRDDATA_H_
#define NRRDDATA_H_

#include "ISignalData.h"

#include <iostream>
#include <string>
#include <vector>
#include <teem/nrrd.h>
#include "linalg.h"

/**
 * \class NrrdData
 * \implements ISignalData
 * \brief Class wrapping teem nrrd for loading nrrd files.
*/

class NrrdData : public ISignalData
{
public:

  /** Constructor */
  NrrdData(ukfPrecisionType sigma_signal, ukfPrecisionType sigma_mask);

  /** Destructor */
  ~NrrdData();

  /** Interpolates the DWI signal at a certain position */
  virtual void Interp3Signal(const vec3_t& pos, ukfVectorType& signal) const;

  /** Interpolates the brain mask at a certain position */
  virtual ukfPrecisionType Interp3ScalarMask(const vec3_t& pos) const;

  /** Gets brain mask value at a certain position */
  virtual ukfPrecisionType ScalarMaskValue(const vec3_t& pos) const;

  /**
   * \brief Get the seed points from the nrrd file
   *
   * Takes care of different seed data types by type casting
   *
   * \param[in]  labels  a vector of labels that define the seed region
   * \param[out] seeds   a vector containing the positions in ijk-space of the seeds
  */
  virtual void GetSeeds(const std::vector<int>& labels, stdVec_t& seeds) const;

  /** returns the gradients of the diffusion image */
  virtual const stdVec_t & gradients() const
  {
    return _gradients;
  }

  /**
   * returns the vector of b-values of the diffusion image<br>
   * Note: Except for cases recorded with multiple b-values it
   *       contains identical values
  */
  virtual const ukfVectorType & GetBValues() const
  {
    return _b_values;
  }

  /**
   * returns the dimension of the signal <br>
   * Note: The actual signal vector will be twice this size
  */
  virtual int GetSignalDimension() const
  {
    return _num_gradients;
  }

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
                        const bool normalizedDWIData, const bool outputNormalizedDWIData);


  virtual bool SetData(Nrrd* data, Nrrd* seed, Nrrd* mask, bool normalizedDWIData);

  /** Returns the dimensions of the signal in each directions as a vector */
  vec3_t dim() const
  {
    return _dim;
  }

private:
  /**
    * Load the signal, called by LoadData
    * \todo Should be a private function of this class, and not implementing ISignalData
  */
  bool LoadSignal(Nrrd* input_nrrd, const bool normalizedDWIData);

  /** The volume dimensions */
  vec3_t _dim;

  int _num_gradients;

  /** gradient directions of the diffusion image */
  stdVec_t _gradients;

  ukfVectorType _b_values;

  /** pointer diffusion data as float */
  float *_data;
  /** pointer to seed data, is casted at runtime */
  void *_seed_data;
  /** seed type is needed for correct casting type */
  int _seed_data_type;
  /** pointer to mask data, is casted at runtime */
  void *_mask_data;
  /** number of bytes of the mask is needed for casting */
  int _mask_num_bytes;

  /** the actual diffusion data in Nrrd type */
  Nrrd *_data_nrrd;
  /** The actual seed data */
  Nrrd *_seed_nrrd;
  /** The actual mask data */
  Nrrd *_mask_nrrd;
};

#endif  // NRRDDATA_H_
