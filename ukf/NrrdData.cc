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
    _data(NULL), _data_nrrd(NULL), _seed_data(NULL), _mask_data(NULL)
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
  ukfPrecisionType       value;

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
      std::cout << "Unsupported data type for seed file!" << std::endl;
      throw;
    }

  return value;
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

void NrrdData::GetSeeds(const std::vector<int>& labels,
                        stdVec_t& seeds) const
{
  if( _seed_data )
    {
    assert(seeds.size() == 0);
    std::vector<int>::const_iterator cit;

    // Go through the volume.
    int nx = _seed_nrrd->axis[2].size;
    int ny = _seed_nrrd->axis[1].size;
    int nz = _seed_nrrd->axis[0].size;
    assert(_seed_data);

    if ( !(nx == _dim[0] && ny == _dim[1] && nz == _dim[2]) )
      {
      itkGenericExceptionMacro(<< "Labelmap ROI volume dimensions DO NOT match DWI dimensions");
      }

    for( int i = 0; i < nx; ++i )
      {
      for( int j = 0; j < ny; ++j )
        {
        for( int k = 0; k < nz; ++k )
          {
          for( cit = labels.begin(); cit != labels.end(); ++cit )
            {
            int value = 0;
            int index = ny * nz * i + nz * j + k;

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
              seeds.push_back(vec3_t(i, j, k) );
              }
            }
          }
        }
      }
    }
  else
    {
    std::cout << "No seed data available." << std::endl;
    }
}

bool NrrdData::SetData(Nrrd* data_nrrd, Nrrd* mask_nrrd, Nrrd* seed_nrrd,
                       bool normalizedDWIData)
  {
  //_data_nrrd = (Nrrd*)data;
  //_seed_data = (Nrrd*)seed;
  //_mask_data = (Nrrd*)mask;

  if( LoadSignal(data_nrrd, normalizedDWIData) )
    {
    return true;
    }

  if (seed_nrrd)
    {
    this->_seed_nrrd = seed_nrrd;
    this->_seed_data = seed_nrrd->data;
    this->_seed_data_type = seed_nrrd->type;
    assert(_seed_data_type == 2 || _seed_data_type == 3 || _seed_data_type == 5);
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

  return false;
  }

bool NrrdData::LoadData(const std::string& data_file,
                        const std::string& seed_file,
                        const std::string& mask_file,
                        const bool normalizedDWIData,
                        const bool outputNormalizedDWIData
                        )
{
  if( _data || _seed_data || _mask_data )
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
  // Load seeds
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

  Nrrd* mask_nrrd = nrrdNew();
  // Load mask
  if( nrrdLoad(mask_nrrd, mask_file.c_str(), NULL) )
    {
    char *err = biffGetDone(NRRD);
    std::cout << "Trouble reading " << mask_file << ": " << err << std::endl;
    free( err );
    return true;
    }

  bool status = SetData(input_nrrd, mask_nrrd, seed_nrrd, normalizedDWIData);
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

  ukfPrecisionType              bValue = ukfZero;

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
  const unsigned int gradientCount = _gradients.size();
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

  _num_gradients = _data_nrrd->axis[0].size;
  assert(_num_gradients == static_cast<int>(gradientCount ) );

  // Get the measurement frame.
  ukfMatrixType measurement_frame(3, 3);
  for( int i = 0; i < 3; ++i )
    {
    for( int j = 0; j < 3; ++j )
      {
//       measurement_frame(i, j) = _data_nrrd->measurementFrame[j][i] ;	//Notice this transpose
      measurement_frame(i, j) = _data_nrrd->measurementFrame[i][j];
      // The Nrrd structure stores the measurement frame vectors as rows in its measurementFrame array
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
    const unsigned int dupIndex = i + gradientCount;
    // Duplicate and reverse direction.
    _gradients.push_back(-_gradients[i]);
    _b_values[dupIndex] =_b_values[i];
    }
  return false;
}
