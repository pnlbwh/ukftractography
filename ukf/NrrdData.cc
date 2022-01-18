/**
 * \file NrrdData.cc
 * \brief implementation of NrrdData.h
*/

#include "NrrdData.h"
#include "ISignalData.h"
#include "dwi_normalize.h"
#include <iostream>
#include <cassert>
#include "itkMacro.h"

NrrdData::NrrdData(ukfPrecisionType sigma_signal, ukfPrecisionType sigma_mask)
  : ISignalData(sigma_signal, sigma_mask),
    _data(NULL), _seed_data(NULL), _stop_data(NULL), _mask_data(NULL), _wm_data(NULL), _gm_data(NULL), _csf_data(NULL), _data_nrrd(NULL)
{

}

NrrdData::~NrrdData()
{
  if( _data_nrrd )
    {
    nrrdNuke(_data_nrrd);
    if( _seed_data )
      {
      nrrdNuke(_seed_nrrd);
      }
    if( _stop_data )
      {
      nrrdNuke(_stop_nrrd);
      }
    if( _wm_data )
      {
      nrrdNuke(_wm_nrrd);
      }
    if( _gm_data )
      {
      nrrdNuke(_gm_nrrd);
      }
    if( _csf_data )
      {
      nrrdNuke(_csf_nrrd);
      }
    if( _mask_data )
      {
      nrrdNuke(_mask_nrrd);
      }
    }
}

void NrrdData::Interp3Signal(const vec3_t& pos,
                             ukfVectorType& signal) const
{
  const int nx = static_cast<const int>(_dim[0]);
  const int ny = static_cast<const int>(_dim[1]);
  const int nz = static_cast<const int>(_dim[2]);

  // If sigmaSignal is not set minimum of voxel size is used for interpolation
  ukfPrecisionType sigma = _sigma_signal;
  if (sigma == 0)
    {
    sigma = std::min(std::min(_voxel[0], _voxel[1]), _voxel[2]);
    }

  ukfPrecisionType w_sum = 1e-16; // this != 0 also doesnt seem to be the problem

  assert(signal.size() == static_cast<unsigned int>(_num_gradients * 2) );
  assert(_data);
  // Is this really necessary?
  for( int i = 0; i < 2 * _num_gradients; ++i )
    {
    signal[i] = ukfZero;
    }

  const int step1 = nz * ny * _num_gradients;
  const int step2 = nz * _num_gradients;
  // for each location
  for( int xx = -1; xx <= 1; ++xx )
    {
    const int x = static_cast<const int>(round(pos[0]) + xx);
    if( x < 0 || nx <= x )
      {
      continue;
      }
    const ukfPrecisionType dx = (x - pos[0]) * _voxel[0];
    const ukfPrecisionType dxx = dx * dx;
    for( int yy = -1; yy <= 1; ++yy )
      {
      const int y = static_cast<const int>(round(pos[1]) + yy);
      if( y < 0 || ny <= y )
        {
        continue;
        }

      const ukfPrecisionType dy = (y - pos[1]) * _voxel[1];
      const ukfPrecisionType dyy = dy * dy;
      for( int zz = -1; zz <= 1; ++zz )
        {
        const int z = static_cast<const int>(round(pos[2]) + zz);
        if( z < 0 || nz <= z )
          {
          continue;
          }
        const ukfPrecisionType dz = (z - pos[2]) * _voxel[2];
        const ukfPrecisionType dzz = dz * dz;

        // gaussian smoothing
        const ukfPrecisionType w = std::exp(-(dxx + dyy + dzz) / sigma);
        // for each gradient direction
        for( int i = 0; i < _num_gradients; ++i )
          {
          // interpolate from all six directions
          signal[i] += w * _data[step1 * x + step2 * y + z * _num_gradients + i];
          }
        // sum of all weights
        w_sum += w;
        }
      }
    }

  // Deleted by Wendy
  // signal shouldn't be halved due to ukfPrecisionType occurance of the gradients
  // reinserted by CB, in order to match the MATLAB code.
  // CB: needs to removed in order to interpolate the signal correctly.
  //w_sum *= 2; // Double each occurance.
  for( int i = 0; i < _num_gradients; ++i )
    {
    signal[i] /= w_sum;

    // Push into second spot.
    signal[i + _num_gradients]  = signal[i];  // Duplicate the signals
    }

}

ukfPrecisionType NrrdData::ScalarMaskValue(const vec3_t& pos) const
{
  const int nx = static_cast<const int>(_dim[0]);
  const int ny = static_cast<const int>(_dim[1]);
  const int nz = static_cast<const int>(_dim[2]);

  unsigned int index;
  ukfPrecisionType value;

  const int x = static_cast<const int>(round(pos[0]));
  const int y = static_cast<const int>(round(pos[1]));
  const int z = static_cast<const int>(round(pos[2]));

  if( (x < 0 || nx <= x) ||
      (y < 0 || ny <= y) ||
      (z < 0 || nz <= z)  )
    {
    return ukfZero;
    }

  index = nz * ny * x + nz * y + z;

  // signed or unsigned doesn't make a difference since masks don't contain any negative values
  switch( _mask_num_bytes )
    {
    case 1:
      {
      value = static_cast<char *>(_mask_data)[index];
      }
      break;
    case 2:
      {
      value = static_cast<short *>(_mask_data)[index];
      }
      break;
    default:
      std::cout << "Unsupported data type for mask file!" << std::endl;
      throw;
    }

  return value;
}

ukfPrecisionType NrrdData::isGMCSFProvided() const
{
  if ( _gm_data || _csf_data || _stop_data)
    {
      return ukfOne;
    }
  else
    {
      return ukfZero;
    }
}

ukfPrecisionType NrrdData::ScalarStopValue(const std::vector<int>& stopLabels, const ukfPrecisionType gm_prob_threshold, const ukfPrecisionType csf_prob_threshold, const vec3_t& pos) const
{
  const int nx = static_cast<const int>(_dim[0]);
  const int ny = static_cast<const int>(_dim[1]);
  const int nz = static_cast<const int>(_dim[2]);

  unsigned int index;
  ukfPrecisionType is_stopping_voxel = ukfZero;

  const int x = static_cast<const int>(round(pos[0]));
  const int y = static_cast<const int>(round(pos[1]));
  const int z = static_cast<const int>(round(pos[2]));

  if( (x < 0 || nx <= x) ||
      (y < 0 || ny <= y) ||
      (z < 0 || nz <= z)  )
    {
    return is_stopping_voxel;
    }

  index = nz * ny * x + nz * y + z;

  if ( _stop_data )
    {

    size_t nx_stop = _stop_nrrd->axis[2].size;
    size_t ny_stop = _stop_nrrd->axis[1].size;
    size_t nz_stop = _stop_nrrd->axis[0].size;
    assert(_stop_data);

    if ( !(nx_stop == _dim[0] && ny_stop == _dim[1] && nz_stop == _dim[2]) )
      {
      itkGenericExceptionMacro(<< "Stopping Labelmap ROI volume dimensions DO NOT match DWI dimensions");
      }

    // signed or unsigned doesn't make a difference since masks don't contain any negative values
    int stop_value = 0;
    switch( _stop_data_type )
      {
      case 0:
        {
        stop_value = 0;
        }
        break;
      case 2:
        {
        stop_value = static_cast<unsigned char *>(_stop_data)[index];
        }
        break;
      case 3:
        {
        stop_value = static_cast<short *>(_stop_data)[index];
        }
        break;
      case 5:
        {
        stop_value = static_cast<int *>(_stop_data)[index];
        }
        break;
      default:
        std::cout << "Unsupported data type for stop file!" << std::endl;
        throw;
      }

    std::vector<int>::const_iterator cit;
    for( cit = stopLabels.begin(); cit != stopLabels.end(); ++cit ) 
      {
        if( *cit == stop_value )
          {
            is_stopping_voxel = ukfOne;
          }
      }

    }
  else // If stopping label map is provided, GM and CSF will not be used.
    {
    if ( _gm_data )
      {
      size_t nx_gm = _gm_nrrd->axis[2].size;
      size_t ny_gm = _gm_nrrd->axis[1].size;
      size_t nz_gm = _gm_nrrd->axis[0].size;
      assert(_gm_data);

      if ( !(nx_gm == _dim[0] && ny_gm == _dim[1] && nz_gm == _dim[2]) )
        {
        itkGenericExceptionMacro(<< "GM segmentation volume dimensions DO NOT match DWI dimensions");
        }

      float gm_prob = 0.0;
      switch( _gm_data_type )
      {
      case 0:
        {
        gm_prob = 0.0;
        }
        break;
      case 9:
        {
        gm_prob = static_cast<float *>(_gm_data)[index];
        }
        break;
      case 10:
        {
        gm_prob = static_cast<float *>(_gm_data)[index];
        }
        break;
      default:
        std::cout << "Unsupported data type for GM segmentation file!" << std::endl;
        throw;
      }

      if (gm_prob > gm_prob_threshold)
        {
          is_stopping_voxel = ukfOne;
        }
      }

    if ( _csf_data )
      {
      size_t nx_csf = _csf_nrrd->axis[2].size;
      size_t ny_csf = _csf_nrrd->axis[1].size;
      size_t nz_csf = _csf_nrrd->axis[0].size;
      assert(_csf_data);

      if ( !(nx_csf == _dim[0] && ny_csf == _dim[1] && nz_csf == _dim[2]) )
        {
        itkGenericExceptionMacro(<< "CSF segmentation volume dimensions DO NOT match DWI dimensions");
        }

      float csf_prob = 0.0;
      switch( _csf_data_type )
      {
      case 0:
        {
        csf_prob = 0.0;
        }
        break;
      case 9:
        {
        csf_prob = static_cast<float *>(_csf_data)[index];
        }
        break;
      case 10:
        {
        csf_prob = static_cast<float *>(_csf_data)[index];
        }
        break;
      default:
        std::cout << "Unsupported data type for GM segmentation file!" << std::endl;
        throw;
      }

      if (csf_prob > csf_prob_threshold)
        {
          is_stopping_voxel = ukfOne;
        }
      }
    }

  return is_stopping_voxel;
}

ukfPrecisionType NrrdData::Interp3ScalarMask(const vec3_t& pos) const
{
  const int nx = static_cast<const int>(_dim[0]);
  const int ny = static_cast<const int>(_dim[1]);
  const int nz = static_cast<const int>(_dim[2]);

  unsigned int index;
  ukfPrecisionType       value;

  ukfPrecisionType w_sum = 1e-16;
  ukfPrecisionType s = ukfZero;

  for( int xx = -1; xx <= 1; xx++ )
    {
    const int x = static_cast<const int>(round(pos[0]) + xx);
    if( x < 0 || nx <= x )
      {
      continue;
      }
    ukfPrecisionType dx = (x - pos[0]) * _voxel[0];
    ukfPrecisionType dxx = dx * dx;
    for( int yy = -1; yy <= 1; yy++ )
      {
      const int y = static_cast<const int>(round(pos[1]) + yy);
      if( y < 0 || ny <= y )
        {
        continue;
        }
      ukfPrecisionType dy = (y - pos[1]) * _voxel[1];
      ukfPrecisionType dyy = dy * dy;
      for( int zz = -1; zz <= 1; zz++ )
        {
        const int z = static_cast<const int>(round(pos[2]) + zz);
        if( z < 0 || nz <= z )
          {
          continue;
          }
        ukfPrecisionType dz = (z - pos[2]) * _voxel[2];
        ukfPrecisionType dzz = dz * dz;

        index = nz * ny * x + nz * y + z;

        // signed or unsigned doesn't make a difference since masks don't contain any negative values
        switch( _mask_num_bytes )
          {
          case 1:
            {
            value = static_cast<char *>(_mask_data)[index];
            }
            break;
          case 2:
            {
            value = static_cast<short *>(_mask_data)[index];
            }
            break;
          default:
            std::cout << "Unsupported data type for seed file!" << std::endl;
            throw;
          }

        ukfPrecisionType w = std::exp(-(dxx + dyy + dzz) / _sigma_mask);
        if( value )
          {
          s += w;
          }

        w_sum += w;
        }
      }
    }

  return s / w_sum;
}

void NrrdData::GetSeeds(const std::vector<int>& labels, const ukfPrecisionType wm_prob_threshold,
                        stdVec_t& seeds) const
{

  const int nx = static_cast<const int>(_dim[0]);
  const int ny = static_cast<const int>(_dim[1]);
  const int nz = static_cast<const int>(_dim[2]);

  if( _seed_data )
    {
      std::cout << "Seeding using Option 3, from a seeding label map." << std::endl;
      size_t nx_seed = _seed_nrrd->axis[2].size;
      size_t ny_seed = _seed_nrrd->axis[1].size;
      size_t nz_seed = _seed_nrrd->axis[0].size;
      assert(_seed_data);

      if ( !(nx_seed == _dim[0] && ny_seed == _dim[1] && nz_seed == _dim[2]) )
        {
        itkGenericExceptionMacro(<< "Seeding Labelmap ROI volume dimensions DO NOT match DWI dimensions");
        }
    }
  else if( _wm_data )
    {
      std::cout << "Seeding using Option 2, from a WM segmentation map." << std::endl;;
      size_t nx_wm = _wm_nrrd->axis[2].size;
      size_t ny_wm = _wm_nrrd->axis[1].size;
      size_t nz_wm = _wm_nrrd->axis[0].size;
      assert(_wm_data);

      if ( !(nx_wm == _dim[0] && ny_wm == _dim[1] && nz_wm == _dim[2]) )
        {
        itkGenericExceptionMacro(<< "WM segmentation volume dimensions DO NOT match DWI dimensions");
        }
    }

  assert(seeds.size() == 0); // Make sure seeds is empty

  for( int i = 0; i < nx; ++i )
    {
    for( int j = 0; j < ny; ++j )
      {
      for( int k = 0; k < nz; ++k )
        {
          int is_seeding_voxel = 0;
          size_t index = ny * nz * i + nz * j + k;

          if( _seed_data )
            {
            std::vector<int>::const_iterator cit;
            for( cit = labels.begin(); cit != labels.end(); ++cit )
              {
              int value = 0;

              switch( _seed_data_type )
                {
                case 2:
                  {
                  value = static_cast<unsigned char *>(_seed_data)[index];
                  }
                  break;
                case 3:
                  {
                  value = static_cast<short *>(_seed_data)[index];
                  }
                  break;
                case 5:
                  {
                  value = static_cast<int *>(_seed_data)[index];
                  }
                  break;
                default:
                  std::cout << "Unsupported data type for seed file!" << std::endl;
                  assert(false);
                }
              if( *cit == value )
                {
                is_seeding_voxel = 1;
                }
              }
            }
          else if ( _wm_data )
            {

            float wm_prob = 0.0;

            switch( _wm_data_type )
              {
              case 0:
                {
                wm_prob = 0.0;
                }
                break;
              case 9:
                {
                wm_prob = static_cast<float *>(_wm_data)[index];
                }
                break;
              case 10:
                {
                wm_prob = static_cast<float *>(_wm_data)[index];
                }
                break;
              default:
                std::cout << "Unsupported data type for WM segmentation file!" << std::endl;
                throw;
              }

            if (wm_prob > wm_prob_threshold) // TODO: may support a user provided threshold.
              {
                is_seeding_voxel = 1;
              }
            }

          if( is_seeding_voxel == 1 )
            {
            seeds.push_back(vec3_t(i, j, k) );
            }
        }
      }
    }
    
  }

bool NrrdData::SetData(Nrrd* data_nrrd, Nrrd* mask_nrrd, 
                       bool normalizedDWIData, 
                       Nrrd* seed_nrrd, Nrrd* stop_nrrd, 
                       Nrrd* wm_nrrd, Nrrd* gm_nrrd, Nrrd* csf_nrrd)
  {
  //_data_nrrd = (Nrrd*)data;
  //_seed_data = (Nrrd*)seed;
  //_mask_data = (Nrrd*)mask;

  if( LoadSignal(data_nrrd, normalizedDWIData) )
    {
    return true;
    }

  if( mask_nrrd->type == 1 || mask_nrrd->type == 2 )
    {
    this->_mask_num_bytes = 1;
    }
  else if( mask_nrrd->type == 3 || mask_nrrd->type == 4 )
    {
    this->_mask_num_bytes = 2;
    }
  else
    {
    std::cout
    << "This implementation only accepts masks of type 'signed char', 'unsigned char', 'short', and 'unsigned short'\n";
    std::cout << "Convert your mask using 'unu convert' and rerun.\n";
    exit(1);
    }

  this->_mask_nrrd = mask_nrrd;
  this->_mask_data = mask_nrrd->data;

  if (seed_nrrd)
    {
    this->_seed_nrrd = seed_nrrd;
    this->_seed_data = seed_nrrd->data;
    this->_seed_data_type = seed_nrrd->type;
    assert(_seed_data_type == 2 || _seed_data_type == 3 || _seed_data_type == 5);
    }

  if (stop_nrrd)
    {
    this->_stop_nrrd = stop_nrrd;
    this->_stop_data = stop_nrrd->data;
    this->_stop_data_type = stop_nrrd->type;
    assert(_stop_data_type == 2 || _stop_data_type == 3 || _stop_data_type == 5);
    }

  if (wm_nrrd)
    {
    this->_wm_nrrd = wm_nrrd;
    this->_wm_data = wm_nrrd->data;
    this->_wm_data_type = wm_nrrd->type;
    assert(_wm_data_type == 9 || _wm_data_type == 10);
    }

  if (gm_nrrd)
    {
    this->_gm_nrrd = gm_nrrd;
    this->_gm_data = gm_nrrd->data;
    this->_gm_data_type = gm_nrrd->type;
    assert(_gm_data_type == 9 || _gm_data_type == 10);
    }

  if (csf_nrrd)
    {
    this->_csf_nrrd = csf_nrrd;
    this->_csf_data = csf_nrrd->data;
    this->_csf_data_type = csf_nrrd->type;
    assert(_csf_data_type == 9 || _csf_data_type == 10);
    }

  return false;
  }

bool NrrdData::LoadData(const std::string& data_file,
                        const std::string& mask_file,
                        const bool normalizedDWIData,
                        const bool outputNormalizedDWIData,
                        const std::string& seed_file,
                        const std::string& stop_file,
                        const std::string& wm_file, 
                        const std::string& gm_file, 
                        const std::string& csf_file)
{
  if( _data || _seed_data || _stop_data || _mask_data || _wm_data || _gm_data || _csf_data )
    {
    std::cout << "There is already some data!" << std::endl;
    return true;
    }

  Nrrd* input_nrrd = nrrdNew();
  if( nrrdLoad(input_nrrd, data_file.c_str(), NULL) )
    {
    char *err = biffGetDone(NRRD);
    std::cout << "Trouble reading " << data_file << ": " << err << std::endl;
    free( err );
    return true;
    }

  Nrrd* mask_nrrd = nrrdNew();
  // Load mask
  if( nrrdLoad(mask_nrrd, mask_file.c_str(), NULL) )
    {
    char *err = biffGetDone(NRRD);
    std::cout << "Trouble reading " << mask_file << ": " << err << std::endl;
    free( err );
    return true;
    }

  if( outputNormalizedDWIData )
    {
    std::string normalizedDataFile = data_file.substr(0, data_file.find_last_of('.') );
    normalizedDataFile.append("_normalized.nrrd");
    std::cout << "Writing normalized signal data to: " << normalizedDataFile << std::endl << std::endl;
    if( nrrdSave(normalizedDataFile.c_str(), _data_nrrd, NULL) )
      {
      std::cout << "Failed while saving the normalized data!" << std::endl;
      char *txt = biffGetDone(NRRD);
      std::cout << txt << std::endl;
      free( txt );
      return true;
      }
    }

  Nrrd* seed_nrrd = nrrdNew();
  // Load seeding map
  if( !seed_file.empty() )
    {
    if( nrrdLoad(seed_nrrd, seed_file.c_str(), NULL) )
      {
      char *err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << seed_file << ": " << err << std::endl;
      free( err );
      return true;
      }
    }

  Nrrd* stop_nrrd = nrrdNew();
  // Load stopping map
  if( !stop_file.empty() )
    {
    if( nrrdLoad(stop_nrrd, stop_file.c_str(), NULL) )
      {
      char *err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << stop_file << ": " << err << std::endl;
      free( err );
      return true;
      }
    }

  Nrrd* wm_nrrd = nrrdNew();
  // Load WM segmentation
  if( !wm_file.empty() )
    {
    if( nrrdLoad(wm_nrrd, wm_file.c_str(), NULL) )
      {
      char *err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << wm_file << ": " << err << std::endl;
      free( err );
      return true;
      }
    }

  Nrrd* gm_nrrd = nrrdNew();
  // Load GM segmentation
  if( !gm_file.empty() )
    {
    if( nrrdLoad(gm_nrrd, gm_file.c_str(), NULL) )
      {
      char *err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << gm_file << ": " << err << std::endl;
      free( err );
      return true;
      }
    }

  Nrrd* csf_nrrd = nrrdNew();
  // Load CSF segmentation
  if( !csf_file.empty() )
    {
    if( nrrdLoad(csf_nrrd, csf_file.c_str(), NULL) )
      {
      char *err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << csf_file << ": " << err << std::endl;
      free( err );
      return true;
      }
    }

  bool status = SetData(input_nrrd, mask_nrrd, normalizedDWIData, seed_nrrd, stop_nrrd, wm_nrrd, gm_nrrd, csf_nrrd);
  return status;
}

bool NrrdData::LoadSignal(Nrrd* input_nrrd, const bool normalizedDWIData)
{
  assert(input_nrrd);

  if( normalizedDWIData )
    {
    this->_data_nrrd = input_nrrd;
    }
  else
    {
    this->_data_nrrd = nrrdNew();
    dwiNormalize(input_nrrd, _data_nrrd);   // Do preprocessing on the data
    }

  // After normalization, the first axis of the nrrd data is the list, namely the gradient axis

  // We might have to extend this later on to support more file formats.
  assert(_data_nrrd->type == 9);

  _data = static_cast<float *>(_data_nrrd->data);

  int len = _data_nrrd->kvpArr->len;
  int size = _data_nrrd->kvpArr->size;

  if( size != 2 )
    {
    size = 2;
    }
  assert(size == 2);

  ukfPrecisionType bValue = ukfZero;

  assert(_gradients.size() == 0);
  // Read key value pairs.
  for( int i = 0; i < len * size; i += size )
    {
    std::string key(_data_nrrd->kvp[i]);

    // std::cout << key << " " << !key.compare("DWMRI_b-value") << std::endl;

    if( !key.compare("DWMRI_b-value") )   // NOTE:compare returns 0 if strings match DWMRI_b-value
      {
      bValue = atof(_data_nrrd->kvp[i + 1]);
      }
    else if( key.length() > 14 &&
             !key.substr(0, 14).compare("DWMRI_gradient") )
      {
      double gx, gy, gz;
      if( 3 != sscanf(_data_nrrd->kvp[i + 1], "%lf %lf %lf", &gx, &gy, &gz) )
        {
        std::cout << "The gradient must have three components!" << std::endl;
        return true;
        }

      // DEBUGING
      //std::cout << "Gradients: " << gx << " " << gy << " " << gz << std::endl;
      _gradients.push_back(vec3_t(gx, gy, gz) );
      }
    else if( !key.compare("modality") )
      {
      assert(!std::string(_data_nrrd->kvp[i + 1]).compare("DWMRI") );
      }
    }
  // if multiple bValues are present the gradient norms are the bValues
  // otherwise the bValues are taken from the header
  // if bValue not in header also take the norm
  // normalizing the gradients
  const size_t gradientCount = _gradients.size();
  _b_values.resize(gradientCount * 2 );
  for( unsigned int i = 0; i < gradientCount; ++i )
    {
    const ukfPrecisionType gradientNorm = _gradients[i].norm();
    const ukfPrecisionType effectiveBvalue = (fabs( gradientNorm - ukfOne ) > 1e-4 ) ? gradientNorm * gradientNorm * bValue: bValue;
    //http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
    //It is after this magnitude rescaling that the nominal bValue (given via "DWMRI_b-value:=bValue") applies.
    _b_values[i] = effectiveBvalue;
    _gradients[i].normalize();
    }

  // Voxel spacing.
  double space_dir[NRRD_SPACE_DIM_MAX];
  double spacing1, spacing2, spacing3;
  nrrdSpacingCalculate(this->_data_nrrd, 1, &spacing1, space_dir);
  nrrdSpacingCalculate(this->_data_nrrd, 2, &spacing2, space_dir);
  nrrdSpacingCalculate(this->_data_nrrd, 3, &spacing3, space_dir);
  _voxel << spacing3, spacing2, spacing1;  // NOTE that the _voxel here is in reverse axis order!

  // make sure something computable is in spacing.
  for(unsigned int i = 0; i < this->_data_nrrd->dim; ++i)
    {
    if(!AIR_EXISTS(space_dir[i]))
      {
      space_dir[i] = 1.0;
      }
    if(!AIR_EXISTS(this->_data_nrrd->spaceOrigin[i]))
      {
      this->_data_nrrd->spaceOrigin[i] = -( (_data_nrrd->axis[i].size / 2) * space_dir[i]);
      }
    }

  // DEBUGING
  // std::cout << "Voxel: " << _voxel[0] << " " << _voxel[1] << " " << _voxel[2] << std::endl;

  // Dimensions
  // NOTICE that the _dim is in reverse axis order!
  _dim << _data_nrrd->axis[3].size, _data_nrrd->axis[2].size,
    _data_nrrd->axis[1].size;

  // std::cout << "dim: " << _dim[0] << " " << _dim[1] << " " << _dim[2] << std::endl;

  _num_gradients = static_cast<int>(_data_nrrd->axis[0].size);
  assert(_num_gradients == static_cast<int>(gradientCount ) );

  // Get the measurement frame.
  ukfMatrixType measurement_frame(3, 3);
  for( int i = 0; i < 3; ++i )
    {
    for( int j = 0; j < 3; ++j )
      {
      // The Nrrd structure stores the measurement frame vector[i] as row[i] in its measurementFrame array
      // So it must be transposed before applying to gradients
      // Note that gradients_ras = measurement_frame * gradients_xyz
      // Other than transposing later, read as (i,j)<--[j][i] in one step
      measurement_frame(i, j) = _data_nrrd->measurementFrame[j][i];

      }
    }

  // Get the ijk->RAS transform matrix
  _i2r.resize(4, 4);
  _i2r.setConstant(ukfZero);

  _i2r(0, 0) = _data_nrrd->axis[1].spaceDirection[0];
  _i2r(1, 0) = _data_nrrd->axis[1].spaceDirection[1];
  _i2r(2, 0) = _data_nrrd->axis[1].spaceDirection[2];
  _i2r(0, 1) = _data_nrrd->axis[2].spaceDirection[0];
  _i2r(1, 1) = _data_nrrd->axis[2].spaceDirection[1];
  _i2r(2, 1) = _data_nrrd->axis[2].spaceDirection[2];
  _i2r(0, 2) = _data_nrrd->axis[3].spaceDirection[0];
  _i2r(1, 2) = _data_nrrd->axis[3].spaceDirection[1];
  _i2r(2, 2) = _data_nrrd->axis[3].spaceDirection[2];
  _i2r(0, 3) = _data_nrrd->spaceOrigin[0];
  _i2r(1, 3) = _data_nrrd->spaceOrigin[1];
  _i2r(2, 3) = _data_nrrd->spaceOrigin[2];
  _i2r(3, 3) = ukfOne;

  // RAS->ijk.
  _r2i = _i2r.inverse();

  // Transform gradients.
  ukfMatrixType R(3, 3);
  R = _i2r.block(0,0,3,3);

  // The gradient should not be affected by voxel size, so factor out the voxel sizes
  // This is equivalent to normalizing the space directions
  const ukfPrecisionType vox_x_inv = ukfOne / _voxel[2];
  const ukfPrecisionType vox_y_inv = ukfOne / _voxel[1];
  const ukfPrecisionType vox_z_inv = ukfOne / _voxel[0];

  R(0, 0) *=  vox_x_inv;
  R(1, 0) *=  vox_x_inv;
  R(2, 0) *=  vox_x_inv;  // R(0,0), R(1,0), R(2,0) is a unit vector, and is just the normalized spacedirection of axis
                          // 1
  R(0, 1) *=  vox_y_inv;
  R(1, 1) *=  vox_y_inv;
  R(2, 1) *=  vox_y_inv;
  R(0, 2) *=  vox_z_inv;
  R(1, 2) *=  vox_z_inv;
  R(2, 2) *=  vox_z_inv;

  ukfMatrixType tmp_mat = R.inverse() * measurement_frame;

  ukfVectorType u(3);
  ukfVectorType u_new(3);
  for( unsigned int i = 0; i < gradientCount; ++i )
    {
    // Transform and normalize.
    u[0] = _gradients[i][0];
    u[1] = _gradients[i][1];
    u[2] = _gradients[i][2];

    u_new = tmp_mat * u;
    const ukfPrecisionType dNorm_inv = ukfOne / u_new.norm();

    // No need to worry about the divison by zero here, since the normalized dwi data has no zero gradient
    _gradients[i][0] = u_new[0] * dNorm_inv;
    _gradients[i][1] = u_new[1] * dNorm_inv;
    _gradients[i][2] = u_new[2] * dNorm_inv;

    }
  // Add reversed gradients
  // This is necessary since the data (signals and gradients) stored in the data files are typically for a half-sphere
  // To get the data for the other half-sphere, simply reverse the gradients and duplicate the signals
  for( unsigned int i = 0; i < gradientCount; ++i )
    {
    const unsigned int dupIndex = static_cast<unsigned int>(i + gradientCount);
    // Duplicate and reverse direction.
    _gradients.push_back(-_gradients[i]);
    _b_values[dupIndex] =_b_values[i];
    }
  return false;
}
