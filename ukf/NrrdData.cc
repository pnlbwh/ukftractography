/**
 * \file NrrdData.cc
 * \brief implementation of NrrdData.h
*/

#include "NrrdData.h"
#include "ISignalData.h"
#include "dwi_normalize.h"
#include <iostream>
#include <cassert>

#include <vnl/vnl_inverse.h>
#include <vnl/vnl_vector.h>

NrrdData::NrrdData(double sigma_signal, double sigma_mask)
  : ISignalData(sigma_signal, sigma_mask), _data(NULL), _seed_data(NULL), _mask_data(NULL)
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

void NrrdData::Interp3Signal(const vec_t& pos,
                             std::vector<double>& signal) const
{
  const int nx = static_cast<const int>(_dim._[0]);
  const int ny = static_cast<const int>(_dim._[1]);
  const int nz = static_cast<const int>(_dim._[2]);

  double w_sum = 1e-16; // this != 0 also doesnt seem to be the problem

  assert(signal.size() == static_cast<size_t>(_num_gradients * 2) );
  assert(_data);
  // Is this really necessary?
  for( int i = 0; i < 2 * _num_gradients; ++i )
    {
    signal[i] = 0.0;
    }

  int step1 = nz * ny * _num_gradients;
  int step2 = nz * _num_gradients;
  // for each location
  for( int xx = -1; xx <= 1; ++xx )
    {
    const int x = static_cast<const int>(round(pos._[0]) + xx);
    if( x < 0 || nx <= x )
      {
      continue;
      }
    double dx = (x - pos._[0]) * _voxel._[0];
    double dxx = dx * dx;
    for( int yy = -1; yy <= 1; ++yy )
      {
      const int y = static_cast<const int>(round(pos._[1]) + yy);
      if( y < 0 || ny <= y )
        {
        continue;
        }

      double dy = (y - pos._[1]) * _voxel._[1];
      double dyy = dy * dy;
      for( int zz = -1; zz <= 1; ++zz )
        {
        const int z = static_cast<const int>(round(pos._[2]) + zz);
        if( z < 0 || nz <= z )
          {
          continue;
          }
        double dz = (z - pos._[2]) * _voxel._[2];
        double dzz = dz * dz;

        // gaussian smoothing
        double w = exp(-(dxx + dyy + dzz) / _sigma_signal);
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
  // signal shouldn't be halved due to double occurance of the gradients
  // reinserted by CB, in order to match the MATLAB code.
  // CB: needs to removed in order to interpolate the signal correctly.
  w_sum *= 2; // Double each occurance.
  for( int i = 0; i < _num_gradients; ++i )
    {
    signal[i] /= w_sum;

    // Push into second spot.
    signal[i + _num_gradients]  = signal[i];  // Duplicate the signals
    }

}

double NrrdData::Interp3ScalarMask(const vec_t& pos) const
{
  const int nx = static_cast<const int>(_dim._[0]);
  const int ny = static_cast<const int>(_dim._[1]);
  const int nz = static_cast<const int>(_dim._[2]);

  unsigned int index;
  double       value;

  double w_sum = 1e-16;
  double s = 0.0;

  for( int xx = -1; xx <= 1; xx++ )
    {
    const int x = static_cast<const int>(round(pos._[0]) + xx);
    if( x < 0 || nx <= x )
      {
      continue;
      }
    double dx = (x - pos._[0]) * _voxel._[0];
    double dxx = dx * dx;
    for( int yy = -1; yy <= 1; yy++ )
      {
      const int y = static_cast<const int>(round(pos._[1]) + yy);
      if( y < 0 || ny <= y )
        {
        continue;
        }
      double dy = (y - pos._[1]) * _voxel._[1];
      double dyy = dy * dy;
      for( int zz = -1; zz <= 1; zz++ )
        {
        const int z = static_cast<const int>(round(pos._[2]) + zz);
        if( z < 0 || nz <= z )
          {
          continue;
          }
        double dz = (z - pos._[2]) * _voxel._[2];
        double dzz = dz * dz;

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
            exit(1);
          }

        double w = exp(-(dxx + dyy + dzz) / _sigma_mask);
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
                        std::vector<vec_t>& seeds) const
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
    assert(nx == _dim._[0] && ny == _dim._[1] && nz == _dim._[2]);
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
              seeds.push_back(make_vec(i, j, k) );
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

  if( LoadSignal(data_file, normalizedDWIData) )
    {
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

  _mask_nrrd = nrrdNew();

  char *err;

  // Load seeds
  if( !seed_file.empty() )
    {
    _seed_nrrd = nrrdNew();
    if( nrrdLoad(_seed_nrrd, seed_file.c_str(), NULL) )
      {
      err = biffGetDone(NRRD);
      std::cout << "Trouble reading " << seed_file << ": " << err << std::endl;
      free( err );
      return true;
      }

    _seed_data_type = _seed_nrrd->type;
    assert(_seed_data_type == 2 || _seed_data_type == 3 || _seed_data_type == 5);
    _seed_data = _seed_nrrd->data;
    }

  // Load mask
  if( nrrdLoad(_mask_nrrd, mask_file.c_str(), NULL) )
    {
    err = biffGetDone(NRRD);
    std::cout << "Trouble reading " << mask_file << ": " << err << std::endl;
    free( err );
    return true;
    }

  if( _mask_nrrd->type == 1 || _mask_nrrd->type == 2 )
    {
    _mask_num_bytes = 1;
    }
  else if( _mask_nrrd->type == 3 || _mask_nrrd->type == 4 )
    {
    _mask_num_bytes = 2;
    }
  else
    {
    std::cout
    << "This implementation only accepts masks of type 'signed char', 'unsigned char', 'short', and 'unsigned short'\n";
    std::cout << "Convert your mask using 'unu convert' and rerun.\n";
    exit(1);
    }

  _mask_data = _mask_nrrd->data;

  return false;
}

bool NrrdData::LoadSignal(const std::string& data_file, const bool normalizedDWIData)
{
  Nrrd *tempNrrd = nrrdNew();

  _data_nrrd = nrrdNew();

  char *err;
  if( nrrdLoad(tempNrrd, data_file.c_str(), NULL) )
    {
    err = biffGetDone(NRRD);
    std::cout << "Trouble reading " << data_file << ": " << err << std::endl;
    free( err );
    return true;
    }

  if( normalizedDWIData )
    {
    nrrdNuke(_data_nrrd);
    _data_nrrd = tempNrrd;
    }
  else
    {
    dwiNormalize(tempNrrd, _data_nrrd);   // Do preprocessing on the data
    nrrdNuke(tempNrrd);
    }
  tempNrrd = NULL;

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

  double              b = 0.0;
  std::vector<double> norm_vec;

  assert(_gradients.size() == 0);
  // Read key value pairs.
  for( int i = 0; i < len * size; i += size )
    {
    std::string key(_data_nrrd->kvp[i]);

    // std::cout << key << " " << !key.compare("DWMRI_b-value") << std::endl;

    if( !key.compare("DWMRI_b-value") )   // NOTE:compare returns 0 if strings match DWMRI_b-value
      {
      b = atof(_data_nrrd->kvp[i + 1]);
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
      // std::cout << "Gradients: " << gx << " " << gy << " " << gz << std::endl;

      // normalizing the gradients

      double norm = sqrt(gx * gx + gy * gy + gz * gz);
      norm_vec.push_back(norm);

      gx /= norm;
      gy /= norm;
      gz /= norm;

      _gradients.push_back(make_vec(gx, gy, gz) );

      }
    else if( !key.compare("modality") )
      {
      assert(!std::string(_data_nrrd->kvp[i + 1]).compare("DWMRI") );
      }
    }

  // if the norm of any gradient norm is >= 1.5 we assume multiple b values
  bool multiple_b_values = false;
  for( unsigned int i = 0; i < _gradients.size(); ++i )
    {
    if( norm_vec[i] > 1.5 )
      {
      multiple_b_values = true;
      break;
      }
    }

  // if multiple b values are present the gradient norms are the b values
  // otherwise the b values are taken from the header
  // if b not in header also take the norm
  _b_values.resize(_gradients.size() );
  for( unsigned int i = 0; i < _gradients.size(); ++i )
    {
    if( multiple_b_values || b == 0.0 )
      {
      _b_values[i] = norm_vec[i];
      }
    else
      {
      _b_values[i] = b;
      }
    }

  // Voxel spacing.
  double space_dir[NRRD_SPACE_DIM_MAX];
  double spacing1, spacing2, spacing3;
  nrrdSpacingCalculate(this->_data_nrrd, 1, &spacing1, space_dir);
  nrrdSpacingCalculate(this->_data_nrrd, 2, &spacing2, space_dir);
  nrrdSpacingCalculate(this->_data_nrrd, 3, &spacing3, space_dir);
  _voxel = make_vec(spacing3, spacing2, spacing1);  // NOTE that the _voxel here is in reverse axis order!

  // DEBUGING
  // std::cout << "Voxel: " << _voxel._[0] << " " << _voxel._[1] << " " << _voxel._[2] << std::endl;

  // Dimensions
  // NOTICE that the _dim is in reverse axis order!
  _dim = make_vec(_data_nrrd->axis[3].size, _data_nrrd->axis[2].size,
                  _data_nrrd->axis[1].size);

  // std::cout << "dim: " << _dim._[0] << " " << _dim._[1] << " " << _dim._[2] << std::endl;

  _num_gradients = _data_nrrd->axis[0].size;
  assert(_num_gradients == static_cast<int>(_gradients.size() ) );

  // Get the measurement frame.
  vnl_matrix<double> measurement_frame(3, 3);
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
  _i2r.set_size(4, 4);
  _i2r.fill(0);

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
  _i2r(3, 3) = 1.0;

  // RAS->ijk.
  _r2i = vnl_inverse(_i2r);

  // Transform gradients.
  vnl_matrix<double> R(3, 3);
  R = _i2r.extract(3, 3);

  // The gradient should not be affected by voxel size, so factor out the voxel sizes
  // This is equivalent to normalizing the space directions
  double vox_x_inv = 1.0 / _voxel._[2];
  double vox_y_inv = 1.0 / _voxel._[1];
  double vox_z_inv = 1.0 / _voxel._[0];

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

  vnl_matrix<double> tmp_mat = vnl_inverse(R) * measurement_frame;

  int                count = _gradients.size();
  vnl_vector<double> u(3), u_new(3);
  for( int i = 0; i < count; ++i )
    {
    // Transform and normalize.
    u[0] = _gradients[i]._[0];
    u[1] = _gradients[i]._[1];
    u[2] = _gradients[i]._[2];

    u_new = tmp_mat * u;
    double dNorm_inv = 1.0 / u_new.magnitude();

    // No need to worry about the divison by zero here, since the normalized dwi data has no zero gradient
    _gradients[i]._[0] = u_new[0] * dNorm_inv;
    _gradients[i]._[1] = u_new[1] * dNorm_inv;
    _gradients[i]._[2] = u_new[2] * dNorm_inv;

    }
  // Add reversed gradients
  // This is necessary since the data (signals and gradients) stored in the data files are typically for a half-sphere
  // To get the data for the other half-sphere, simply reverse the gradients and duplicate the signals
  for( int i = 0; i < count; ++i )
    {
    // Duplicate and reverse direction.
    _gradients.push_back(-_gradients[i]);
    _b_values.push_back(_b_values[i]);
    }
  return false;
}
