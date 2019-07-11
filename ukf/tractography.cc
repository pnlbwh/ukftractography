/**
 * \file tractography.cc
 * \brief implementation of tractography.h
*/

#include <itkMacro.h> // needed for ITK_VERSION_MAJOR
#if ITK_VERSION_MAJOR >= 5
#include "itkMultiThreaderBase.h"
#include "itkPlatformMultiThreader.h"
#else

#include "itkMultiThreader.h"
#endif

// System includes
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>

// VTK includes
#include "vtkNew.h"
#include "vtkPolyData.h"

// UKF includes
#include "tractography.h"
#include "filter_model.h"
#include "ISignalData.h"
#include "NrrdData.h"
#include "utilities.h"
#include "vtk_writer.h"
#include "thread.h"
#include "math_utilities.h"

// filters
#include "filter_Full1T.h"
#include "filter_Full1T_FW.h"
#include "filter_Full2T.h"
#include "filter_Full2T_FW.h"
#include "filter_Full3T.h"
#include "filter_NODDI1F.h"
#include "filter_NODDI2F.h"
#include "filter_Simple1T.h"
#include "filter_Simple1T_FW.h"
#include "filter_Simple2T.h"
#include "filter_Simple2T_FW.h"
#include "filter_Simple3T.h"

// TODO implement this switch
//#include "config.h"

Tractography::Tractography(UKFSettings s) :

    // begin initializer list

    _ukf(0, NULL),



    _output_file(s.output_file),
    _output_file_with_second_tensor(s.output_file_with_second_tensor),

    _record_fa    (s.record_fa),
    _record_nmse  (s.record_nmse),
    _record_trace (s.record_trace),
    _record_state (s.record_state),
    _record_cov   (s.record_cov),
    _record_free_water(s.record_free_water),
    _record_Vic   (s.record_Vic),
    _record_kappa (s.record_kappa),
    _record_Viso  (s.record_Viso),
    _record_tensors(s.record_tensors),
    _transform_position(s.transform_position),
    _store_glyphs (s.store_glyphs),
    _branches_only(s.branches_only),

    _p0(s.p0),
    _sigma_signal(s.sigma_signal),
    _sigma_mask(s.sigma_mask),
    _min_radius(s.min_radius),
    _max_length(static_cast<int>(std::ceil(s.maxHalfFiberLength / s.stepLength) ) ),
    _full_brain(false),
    _noddi(s.noddi),
    _fa_min(s.fa_min), _mean_signal_min(s.mean_signal_min),
    _seeding_threshold(s.seeding_threshold),
    _num_tensors(s.num_tensors),
    _seeds_per_voxel(s.seeds_per_voxel),
    _cos_theta_min(std::cos(DegToRad(s.min_branching_angle))),
    _cos_theta_max(std::cos(DegToRad(s.max_branching_angle))),
    _is_full_model(s.is_full_model),
    _free_water(s.free_water),
    _stepLength(s.stepLength),
    _steps_per_record(s.recordLength/s.stepLength),
    _labels(s.labels),

    Qm(s.Qm),
    Ql(s.Ql),
    Qw(s.Qw),
    Qkappa(s.Qkappa),
    Qvic(s.Qvic),
    Rs(s.Rs),

    _writeBinary(true),
    _writeCompressed(true),

    _num_threads(s.num_threads),
    _outputPolyData(NULL),

    _filter_model_type(Tractography::_1T),
    _model(NULL),
    debug(false)
    // end initializer list
{
}

Tractography::~Tractography()
{
  if( this->_signal_data )
    {
    delete this->_signal_data;
    }
  if (this->_model)
  {
    delete this->_model; // TODO smartpointer
  }
}

void Tractography::UpdateFilterModelType()
{
  if (!this->_signal_data)
    {
    return;
    }

  if (this->_model) {
    delete this->_model; // TODO smartpointer
  }

  this->_filter_model_type = Tractography::_1T;
  bool simpleTensorModel = !(this->_is_full_model);

  if (_noddi){
    if (_num_tensors == 1){
      std::cout << "Using NODDI 1-Fiber model." << std::endl;
      this->_filter_model_type = Tractography::_1T_FW; // same vtk writer can be used
    } else if (_num_tensors == 2){
      std::cout << "Using NODDI 2-Fiber model." << std::endl;
      this->_filter_model_type = Tractography::_2T_FW; // same vtk writer can be used
    }
  } else if (_num_tensors == 1) {
    if (simpleTensorModel && !_free_water) {
      std::cout << "Using 1-tensor simple model." << std::endl;
      this->_filter_model_type = Tractography::_1T;
    } else if (simpleTensorModel && _free_water) {
      std::cout << "Using 1-tensor simple model with free water estimation." << std::endl;
      this->_filter_model_type = Tractography::_1T_FW;
    } else if (!simpleTensorModel && !_free_water) {
      std::cout << "Using 1-tensor full model." << std::endl;
      this->_filter_model_type = Tractography::_1T_FULL;
    } else if (!simpleTensorModel && _free_water) {
      std::cout << "Using 1-tensor full model with free water estimation." << std::endl;
      this->_filter_model_type = Tractography::_1T_FW_FULL;
    }
  } else if (_num_tensors == 2) {
    if (simpleTensorModel && !_free_water) {
      std::cout << "Using 2-tensor simple model." << std::endl;
      this->_filter_model_type = Tractography::_2T;
    } else if (simpleTensorModel && _free_water) {
      std::cout << "Using 2-tensor simple model with free water estimation." << std::endl;
      this->_filter_model_type = Tractography::_2T_FW;
    } else if (!simpleTensorModel && !_free_water) {
      std::cout << "Using 2-tensor full model." << std::endl;
      this->_filter_model_type = Tractography::_2T_FULL;
    } else if (!simpleTensorModel && _free_water) {
      std::cout << "Using 2-tensor full model with free water estimation." << std::endl;
      this->_filter_model_type = Tractography::_2T_FW_FULL;
    }
  } else if (_num_tensors == 3) {
    if (simpleTensorModel) {
      std::cout << "Using 3-tensor simple model." << std::endl;
      this->_filter_model_type = Tractography::_3T;
    } else {
      std::cout << "Using 3-tensor full model." << std::endl;
      this->_filter_model_type = Tractography::_3T_FULL;
    }
  }

  if( (_cos_theta_max != ukfOne) && (_cos_theta_max <= _cos_theta_min) )
    {
    std::cout << "Maximum branching angle must be greater than " << RadToDeg(_cos_theta_min) << " degrees." << std::endl;
    throw;
    }

  if( _num_tensors < 1 || _num_tensors > 3 )
    {
    std::cout << "Only one, two or three tensors are supported." << std::endl;
    throw;
    }

  if (_noddi) {
    if (_record_fa || _record_trace || _record_free_water || _record_tensors) {
        std::cout << "recordFA, recordTrace, record_free_water, recordTensors parameters can only be used with tensor models\n";
        throw;
    }
  }

  // Double check branching.
  _is_branching = _num_tensors > 1 && _cos_theta_max < ukfOne;  // The branching is enabled when the maximum branching
                                                             // angle is not 0
  std::cout << "Branching " << (_is_branching ? "enabled" : "disabled") << std::endl << std::endl;
  if( !_is_branching )
    {
    _branches_only = false;
    }

  _nPosFreeWater = -1; // not used for normal case
  // for free water case used in the Record function to know where the fw is in the state
  if( _num_tensors == 1 )   // 1 TENSOR CASE /////////////////////////////////////
    {
    if( !_is_full_model && _free_water ) // simple model with free water
      {
      _nPosFreeWater = 5;
      }
    else if( _is_full_model && _free_water ) // full model with free water
      {
      _nPosFreeWater = 6;
      }
    else if( _noddi ) // Viso is recorded
      {
      _nPosFreeWater = 5;
      }
    }
  else if( _num_tensors == 2 )    // 2 TENSOR CASE ////////////////////////////////
    {
    if( !_is_full_model && _free_water ) // simple model with free water
      {
      _nPosFreeWater = 10;
      }
    else if( _is_full_model && _free_water ) // full model with free water
      {
      _nPosFreeWater = 12;
      }
    else if( _noddi ) // Viso is recorded
      {
      _nPosFreeWater = 10;
      }
    }

  // set up tensor weights
  this->weights_on_tensors.resize(_num_tensors);
  for (int i = 0; i < _num_tensors; i++)
    {
    this->weights_on_tensors[i]=(1.0 / _num_tensors) ;
    }

  ukfPrecisionType weight_accumu = 0 ;
  for (int i = 0; i < _num_tensors; i++)
    {
    weight_accumu += weights_on_tensors[i] ;
    }
  if (std::abs(weight_accumu - 1.0) > 0.000001)
    {
    std::cout << "The weights on different tensors must add up to 1!" << std::endl << std::endl ;
    throw;
    }
  else
    {
    weights_on_tensors.norm(); // Normalize for all to add up to 1.
    }


  // TODO refactor this NODDI switch
  if (this->_filter_model_type == _1T_FW && this->_noddi && this->_num_tensors == 1) {
    _model = new NODDI1F(Qm, Qkappa, Qvic, Rs, this->weights_on_tensors, this->_noddi);
  }
  else if (this->_filter_model_type == _2T_FW && this->_noddi && this->_num_tensors == 2) {
    _model = new NODDI2F(Qm, Qkappa, Qvic, Rs, this->weights_on_tensors, this->_noddi);
  }
  else if (this->_filter_model_type == _1T) {
    _model = new Simple1T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _1T_FW) {
    _model = new Simple1T_FW(Qm, Ql, Qw, Rs, this->weights_on_tensors, this->_free_water, D_ISO);
  }
  else if (this->_filter_model_type == _1T_FULL) {
    _model = new Full1T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _1T_FW_FULL) {
    _model = new Full1T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _2T) {
    _model = new Simple2T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _2T_FW) {
    _model = new Simple2T_FW(Qm, Ql, Qw, Rs, this->weights_on_tensors, this->_free_water, D_ISO);
  }
  else if (this->_filter_model_type == _2T_FULL) {
    _model = new Full2T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _2T_FW_FULL) {
    _model = new Full2T_FW(Qm, Ql, Qw, Rs, this->weights_on_tensors, this->_free_water, D_ISO);
  }
  else if (this->_filter_model_type == _3T) {
    _model = new Simple3T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else if (this->_filter_model_type == _3T_FULL) {
    _model = new Full3T(Qm, Ql, Rs, this->weights_on_tensors, this->_free_water);
  }
  else {
    std::cerr << "Unknown filter type!" << std::endl;
    assert(false);
  }

  _model->set_signal_data(_signal_data);
  _model->set_signal_dim(_signal_data->GetSignalDimension() * 2);
}

bool Tractography::SetData(void* data, void* mask, void* seed,
                           bool normalizedDWIData)
{
  if (!data || !mask)
    {
    std::cout << "Invalid input Nrrd pointers!" << std::endl;
    return true;
    }

  if (!seed)
    {
    _full_brain = true;
    }

  _signal_data = new NrrdData(_sigma_signal, _sigma_mask);
  _signal_data->SetData((Nrrd*)data, (Nrrd*)mask, (Nrrd*)seed, normalizedDWIData);

  return false;
}

bool Tractography::LoadFiles(const std::string& data_file,
                             const std::string& seed_file,
                             const std::string& mask_file,
                             const bool normalized_DWI_data,
                             const bool output_normalized_DWI_data
                             )
{
  _signal_data = new NrrdData(_sigma_signal, _sigma_mask);

  if( seed_file.empty() )
    {
    _full_brain = true;
    }

  if( _signal_data->LoadData(data_file, seed_file, mask_file, normalized_DWI_data, output_normalized_DWI_data) )
    {
    std::cout << "ISignalData could not be loaded" << std::endl;
    delete _signal_data;
    _signal_data = NULL;
    return true;
    }
  return false;
}

void Tractography::Init(std::vector<SeedPointInfo>& seed_infos)
{
  if (! (_signal_data)) {
    std::cout << "No signal data!";
    throw;
  }

  int signal_dim = _signal_data->GetSignalDimension();

  stdVec_t seeds;
  if (! (_labels.size() > 0))
    {
    std::cout << "No label data!";
    throw;
    }

  if(!_ext_seeds.empty())
    {
    seeds = _ext_seeds;
    }
  else if(!_full_brain)
    {
    _signal_data->GetSeeds(_labels, seeds);
    }
  else
    {
    // Iterate through all brain voxels and take those as seeds voxels.
    const vec3_t dim = _signal_data->dim();
    for( int x = 0; x < dim[0]; ++x )
      {
      for( int y = 0; y < dim[1]; ++y )
        {
        for( int z = 0; z < dim[2]; ++z )
          {
          vec3_t pos(x,y,z); //  = make_vec(x, y, z);
          if(_signal_data->ScalarMaskValue(pos) > 0)
            {
            seeds.push_back(pos);
            }
          }
        }
      }
    }

  if (! (seeds.size() > 0)) {
	std::cout << "No matching label ROI seeds found! Please verify label selection.";
    throw;
  }

  // Determinism.
  srand(0);

  // Create random offsets from the seed voxel.
  stdVec_t rand_dirs;

  if( seeds.size() == 1 && _seeds_per_voxel == 1 )   // if there is only one seed don't use offset so fibers can be
                                                     // compared
    {
    rand_dirs.push_back(vec3_t(0,0,0) /* make_vec(0, 0, 0) */);   // in the test cases.
    }
  else
    {
    for( int i = 0; i < _seeds_per_voxel; ++i )
      {
      vec3_t dir(static_cast<ukfPrecisionType>( (rand() % 10001) - 5000),
                static_cast<ukfPrecisionType>( (rand() % 10001) - 5000),
                static_cast<ukfPrecisionType>( (rand() % 10001) - 5000) );

      // CB: those directions are to compare against the matlab output
      // dir._[2] = 0.439598093988175;
      // dir._[1] = 0.236539281163321;
      // dir._[0] = 0.028331682419209;

      dir = dir.normalized();
      dir *= ukfHalf;

      rand_dirs.push_back(dir);
      }
    }
  // Calculate all starting points.
  stdVec_t      starting_points;
  stdEigVec_t   signal_values;
  ukfVectorType signal(signal_dim * 2);

  int num_less_than_zero = 0;
  int num_invalid = 0;
  int num_mean_signal_too_low = 0;

  int tmp_counter = 1;
  for( stdVec_t::const_iterator cit = seeds.begin(); cit != seeds.end(); ++cit )
    {
    for( stdVec_t::iterator jt = rand_dirs.begin(); jt != rand_dirs.end(); ++jt )
      {
      vec3_t point = *cit + *jt;

      _signal_data->Interp3Signal(point, signal); // here and in every step
      tmp_counter++;

      // DEBUG
      // std::cout << "point: " << point._[0] << " " << point._[1] << " " << point._[2] << std::endl;

      // Filter out all starting points that have negative signal values (due to
      // noise) or that otherwise have invalid signal values.
      bool keep = true;
      // We only scan half of the signal values since the second half is simply
      // a copy of the first half.
      for( int k = 0; k < signal_dim; ++k )
        {
        if( signal[k] < 0 )
          {
          keep = false;
          ++num_less_than_zero;
          break;
          }

        if( std::isnan(signal[k]) || std::isinf(signal[k]) )
          {
          keep = false;
          ++num_invalid;
          break;
          }

        // If we do full brain tractography we only want seed voxels where the
        // GA is bigger than 0.18.
        ukfMatrixType signal_tmp(signal_dim * 2, 1);
        signal_tmp.col(0) = signal;
        if( _full_brain && s2adc(signal_tmp) < _seeding_threshold )
          {
          keep = false;
          ++num_mean_signal_too_low;
          break;
          }
        }

      // If all the criteria is met we keep that point and the signal data.
      if( keep )
        {
        signal_values.push_back(signal);
        starting_points.push_back(point);
        }
      }
    }
  stdEigVec_t starting_params(starting_points.size() );

  UnpackTensor(_signal_data->GetBValues(), _signal_data->gradients(),
               signal_values, starting_params);

  // If we work with the simple model we have to change the second and third
  // eigenvalues: l2 = l3 = (l2 + l3) / 2.
  if( !_is_full_model )   // i.e. simple model
    {
    for( size_t i = 0; i < starting_params.size(); ++i )
      {
      starting_params[i][7] = starting_params[i][8] = (starting_params[i][7] + starting_params[i][8]) / 2.0;
      // two minor eigenvalues are treated equal in simplemodel
      }
    }

  // Pack information for each seed point.
  int fa_too_low = 0;
  for( size_t i = 0; i < starting_points.size(); ++i )
    {
    const ukfVectorType& param = starting_params[i];

    assert(param.size() == 9);

    // Filter out seeds whose FA is too low.
    ukfPrecisionType fa = l2fa(param[6], param[7], param[8]);
    ukfPrecisionType trace = param[6] + param[7] + param[8];
    ukfPrecisionType fa2 = -1;
    ukfPrecisionType trace2 = -1;

    if( _num_tensors >= 2 )
      {
      fa2 = fa;
      trace2 = trace;
      }

    if( fa <= _seeding_threshold )
      {
      ++fa_too_low;
      continue;
      }

    // Create seed info for both directions;
    SeedPointInfo info;
    stdVecState tmp_info_state;
    stdVecState tmp_info_inv_state;
    SeedPointInfo info_inv;

    info.point = starting_points[i];
    info.start_dir << param[0], param[1], param[2];
    info.fa = fa;
    info.fa2 = fa2;
    info.trace = trace;
    info.trace2 = trace2;
    info_inv.point = starting_points[i];
    info_inv.start_dir << -param[0], -param[1], -param[2];
    info_inv.fa = fa;
    info_inv.fa2 = fa2;
    info_inv.trace = trace;
    info_inv.trace2 = trace2;

    if( _is_full_model )
      {
      tmp_info_state.resize(6);
      tmp_info_inv_state.resize(6);
      tmp_info_state[0] = param[3];     // Theta
      tmp_info_state[1] = param[4];     // Phi
      tmp_info_state[2] = param[5];     // Psi
      tmp_info_state[5] = param[8];     // l3
      tmp_info_inv_state[0] = param[3]; // Theta
      tmp_info_inv_state[1] = param[4]; // Phi
      // Switch psi angle.
      tmp_info_inv_state[2] = (param[5] < ukfZero ? param[5] + UKF_PI : param[5] - UKF_PI);
      tmp_info_inv_state[5] = param[8]; // l3

      }
    else     // i.e. simple model
      { // Starting direction.
      tmp_info_state.resize(5);
      tmp_info_inv_state.resize(5);
      tmp_info_state[0] = info.start_dir[0];
      tmp_info_state[1] = info.start_dir[1];
      tmp_info_state[2] = info.start_dir[2];
      tmp_info_inv_state[0] = info_inv.start_dir[0];
      tmp_info_inv_state[1] = info_inv.start_dir[1];
      tmp_info_inv_state[2] = info_inv.start_dir[2];
      }

    ukfPrecisionType Viso;
    if(_noddi)
    {
      ukfPrecisionType minnmse,nmse;
      State state(_model->state_dim() );
      minnmse = 99999;
      _signal_data->Interp3Signal(info.point, signal); // here and in every step

      ukfPrecisionType dPar, dIso;
      int n = 5;
      dPar = 0.0000000017;
      dIso = 0.000000003;
      createProtocol(_signal_data->GetBValues(), _gradientStrength, _pulseSeparation);
      double kappas[5] = {0.5, 1, 2, 4, 8};
      double Vic[5] = {0, 0.25, 0.5, 0.75, 1};
      double Visoperm[5] = {0, 0.25, 0.5, 0.75, 1};
      for( int a = 0; a < n; ++a )
      for( int b = 0; b < n; ++b )
      for( int c = 0; c < n; ++c )
      {
        ukfVectorType Eec, Eic, Eiso, E;
        ExtraCelluarModel(dPar, Vic[b], kappas[a], _gradientStrength,
                          _pulseSeparation, _signal_data->gradients(), info.start_dir, Eec);
        IntraCelluarModel(dPar, kappas[a], _gradientStrength, _pulseSeparation,
                          _signal_data->gradients(), info.start_dir, Eic);
        IsoModel(dIso, _gradientStrength, _pulseSeparation, Eiso);

        E.resize(Eic.size());
        E = Visoperm[c]*Eiso + (1 - Visoperm[c])*(Vic[b]*Eic + (1-Vic[b])*Eec);

        nmse = (E-signal).squaredNorm() / signal.squaredNorm();

        if(nmse < minnmse)
        {
          minnmse = nmse;
          Viso = Visoperm[c];
          tmp_info_state[3] = Vic[b];     // Vic
          tmp_info_state[4] = kappas[a];     // Kappa
          tmp_info_inv_state[3] = Vic[b]; // Vic
          tmp_info_inv_state[4] = kappas[a]; // Kappa
       }
      }
      // xstd::cout <<"nmse of initialization "<< minnmse << "\n";
      if(_num_tensors > 1)
      {
        tmp_info_state.resize(10);
        tmp_info_inv_state.resize(10);
        tmp_info_state[5] = info.start_dir[0];
        tmp_info_state[6] = info.start_dir[1];
        tmp_info_state[7] = info.start_dir[2];
        tmp_info_inv_state[5] = info_inv.start_dir[0];
        tmp_info_inv_state[6] = info_inv.start_dir[1];
        tmp_info_inv_state[7] = info_inv.start_dir[2];
        tmp_info_state[8] = tmp_info_state[3];     // Vic
        tmp_info_state[9] = tmp_info_state[4];     // Kappa
        tmp_info_inv_state[8] = tmp_info_inv_state[3]; // Vic
        tmp_info_inv_state[9] = tmp_info_inv_state[4]; //kappa
      }
    }
    else
    {
      tmp_info_state[3] = param[6];     // l1
      tmp_info_state[4] = param[7];     // l2
      tmp_info_inv_state[3] = param[6]; // l1
      tmp_info_inv_state[4] = param[7]; // l2
    }
    // Duplicate/tripple states if we have several tensors.
    if( _num_tensors > 1 && !_noddi)
      {
      size_t size = tmp_info_state.size();
      for( size_t j = 0; j < size; ++j )
        {
        tmp_info_state.push_back(tmp_info_state[j]);
        tmp_info_inv_state.push_back(tmp_info_state[j]);
        }
      if( _num_tensors > 2 )
        {
        for( size_t j = 0; j < size; ++j )
          {
          tmp_info_state.push_back(tmp_info_state[j]);
          tmp_info_inv_state.push_back(tmp_info_state[j]);
          }
        }
      }
    if(_noddi)
    {
      tmp_info_state.push_back(Viso);
      tmp_info_inv_state.push_back(Viso); // add the weight to the state (well was sich rhymt das stiimt)
    }
    else
    {
      if( _free_water )
        {
        tmp_info_state.push_back(1);
        tmp_info_inv_state.push_back(1); // add the weight to the state (well was sich rhymt das stiimt)
        }
    }

    const int state_dim = tmp_info_state.size();

    info.covariance.resize(state_dim, state_dim);
    info_inv.covariance.resize(state_dim, state_dim);

    // make sure covariances are really empty
    info.covariance.setConstant(ukfZero);
    info_inv.covariance.setConstant(ukfZero);
    for( int local_i = 0; local_i < state_dim; ++local_i )
      {
      info.covariance(local_i, local_i) = _p0;
      info_inv.covariance(local_i, local_i) = _p0;
      }
    info.state = ConvertVector<stdVecState,State>(tmp_info_state);
    info_inv.state = ConvertVector<stdVecState,State>(tmp_info_inv_state);
    seed_infos.push_back(info);
    seed_infos.push_back(info_inv);   // NOTE that the seed in reverse direction is put directly after the seed in
                                      // original direction
    }
}

bool Tractography::Run()
{
  assert(_signal_data);   // The _signal_data is initialized in Tractography::LoadFiles(),
  // Thus Run() must be invoked after LoadFiles()
  // Initialize and prepare seeds.

  std::vector<SeedPointInfo>            primary_seed_infos;
  std::vector<SeedPointInfo>            branch_seed_infos;       // The info of branching seeds
  std::vector<BranchingSeedAffiliation> branch_seed_affiliation; // Which fiber originated from the main seeds is this
                                                                 // branch attached

  Init(primary_seed_infos);
  if (primary_seed_infos.size() < 1)
    {
    std::cerr << "No valid seed points available!" << std::endl;
    return false;
    }

  const int num_of_threads = std::min(_num_threads, static_cast<int>(primary_seed_infos.size() ) );
  // const int num_of_threads = 8;

  assert(num_of_threads > 0);
  _ukf.reserve(num_of_threads); //Allocate, but do not assign
  for( int i = 0; i < num_of_threads; i++ )
    {
    _ukf.push_back(new UnscentedKalmanFilter(_model) );   // Create one Kalman filter for each thread
    }

  std::vector<UKFFiber> raw_primary;

    {
    if (this->debug) std::cout << "Tracing " << primary_seed_infos.size() << " primary fibers:" << std::endl;

    raw_primary.resize(primary_seed_infos.size() );

    WorkDistribution work_distribution = GenerateWorkDistribution(num_of_threads,
                                                                  static_cast<int>(primary_seed_infos.size() ) );
#if ITK_VERSION_MAJOR >= 5
    itk::PlatformMultiThreader::Pointer threader = itk::PlatformMultiThreader::New();
    threader->SetNumberOfWorkUnits( num_of_threads );
    std::vector<std::thread> vectorOfThreads;
    vectorOfThreads.reserve(num_of_threads);
#else
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( num_of_threads );
#endif
    thread_struct str;
    str.tractography_ = this;
    str.seed_infos_ = &primary_seed_infos;
    str.work_distribution = &work_distribution;
    str.branching_ = _is_branching;
    str.num_tensors_ = _num_tensors;
    str.output_fiber_group_ = &raw_primary;
    str.branching_seed_info_vec = new std::vector<std::vector<SeedPointInfo> >(num_of_threads);
    str.branching_seed_affiliation_vec = new std::vector<std::vector<BranchingSeedAffiliation> >(num_of_threads);
    for( int i = 0; i < num_of_threads; i++ )
      {
#if ITK_VERSION_MAJOR >= 5
      vectorOfThreads.push_back( std::thread(ThreadCallback, i, &str ) );
#else
      threader->SetMultipleMethod(i, ThreadCallback, &str);
#endif
      }
#if ITK_VERSION_MAJOR < 5
    threader->SetGlobalDefaultNumberOfThreads(num_of_threads);
#else
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(num_of_threads);
#endif

#if ITK_VERSION_MAJOR >= 5
    for(auto & li : vectorOfThreads )
      {
      if(li.joinable())
        {
        li.join();
        }
      }
#else
    threader->MultipleMethodExecute();
#endif

    // Unpack the branch seeds and their affiliation
    int num_branch_seeds = 0;
    for( int i = 0; i < num_of_threads; i++ )
      {
      num_branch_seeds += static_cast<int>(str.branching_seed_info_vec->at(i).size() );
      }

    if (this->debug) std::cout << "branch_seeds size: " << num_branch_seeds << std::endl;
    branch_seed_infos.resize(num_branch_seeds);
    branch_seed_affiliation.resize(num_branch_seeds);

    int counter = 0;
    for( int i = 0; i < num_of_threads; i++ )
      {
      for( size_t j = 0; j < str.branching_seed_info_vec->at(i).size(); j++ )
        {
        branch_seed_infos[counter] = str.branching_seed_info_vec->at(i).at(j);
        branch_seed_affiliation[counter] = str.branching_seed_affiliation_vec->at(i).at(j);
        counter++;
        }
      }
    }

  std::vector<UKFFiber> raw_branch;
  if( _is_branching )
    {
    assert(_num_tensors == 2 || _num_tensors == 3);
    if (this->debug) std::cout << "Tracing " << branch_seed_infos.size() << " branches:" << std::endl;

    raw_branch.resize(branch_seed_infos.size() );

    WorkDistribution work_distribution =
      GenerateWorkDistribution(num_of_threads, static_cast<int>(branch_seed_infos.size() ) );

    thread_struct str;
    str.tractography_ = this;
    str.work_distribution = &work_distribution;
    str.seed_infos_ = &branch_seed_infos;
    str.branching_ = false;
    str.num_tensors_ = _num_tensors;
    str.output_fiber_group_ = &raw_branch;
    str.branching_seed_info_vec = new std::vector<std::vector<SeedPointInfo> >(num_of_threads);
    str.branching_seed_affiliation_vec = new std::vector<std::vector<BranchingSeedAffiliation> >(num_of_threads);

#if ITK_VERSION_MAJOR >= 5
    itk::PlatformMultiThreader::Pointer threader = itk::PlatformMultiThreader::New();
    threader->SetNumberOfWorkUnits( num_of_threads );
    std::vector<std::thread> vectorOfThreads;
    vectorOfThreads.reserve(num_of_threads);
#else
    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( num_of_threads );
#endif


    for( int i = 0; i < num_of_threads; i++ )
      {
#if ITK_VERSION_MAJOR >= 5
      vectorOfThreads.push_back( std::thread( ThreadCallback, i, &str ) );
#else
      threader->SetMultipleMethod(i, ThreadCallback, &str);
#endif
      }
#if ITK_VERSION_MAJOR >= 5
    for(auto & li : vectorOfThreads )
    {
      if(li.joinable())
      {
      li.join();
      }
    }
#else
    threader->MultipleMethodExecute();
#endif
    }

  std::vector<UKFFiber> fibers;
  PostProcessFibers(raw_primary, raw_branch, branch_seed_affiliation, _branches_only, fibers);

  if (this->debug) std::cout << "fiber size after PostProcessFibers: " << fibers.size() << std::endl;

  if ( fibers.size()  == 0 )
  {
    std::cout << "No fibers! Returning." << fibers.size() << std::endl;
    return EXIT_FAILURE;
  }

  // Write the fiber data to the output vtk file.
  VtkWriter writer(_signal_data, this->_filter_model_type, _record_tensors);
  writer.set_transform_position(_transform_position);

  int writeStatus = EXIT_SUCCESS;
  if (this->_outputPolyData != NULL)
  // TODO refactor this is bad control flow
    {
    writer.PopulateFibersAndTensors(this->_outputPolyData, fibers);
    this->_outputPolyData->Modified();
    }
  else
    {
    vtkNew<vtkPolyData> pd;
    this->SetOutputPolyData(pd.GetPointer());
    writer.PopulateFibersAndTensors(this->_outputPolyData, fibers);
    this->_outputPolyData->Modified();

    // possibly write binary VTK file.
    writer.SetWriteBinary(this->_writeBinary);
    writer.SetWriteCompressed(this->_writeCompressed);

    writeStatus = writer.Write(_output_file, _output_file_with_second_tensor,
                                         fibers, _record_state, _store_glyphs, _noddi);

    // TODO refactor!
    this->SetOutputPolyData(NULL);
    }

  // Clear up the kalman filters
  for( size_t i = 0; i < _ukf.size(); i++ )
    {
    delete _ukf[i];
    }
  _ukf.clear();
  return writeStatus;
}

// FIXME: not clear why gradientStrength and pulseSeparation are passed as arguments when
//        they are already class members.
void Tractography::createProtocol(const ukfVectorType& _b_values,
                    ukfVectorType& l_gradientStrength, ukfVectorType& l_pulseSeparation)
{
  std::vector<double> Bunique, tmpG;
  ukfPrecisionType Bmax = 0;
  ukfPrecisionType tmp, Gmax, GAMMA;

  l_gradientStrength.resize(_b_values.size());
  l_pulseSeparation.resize(_b_values.size());

  // set maximum G = 40 mT/m
  Gmax = 0.04;
  GAMMA = 267598700;

  for(int i = 0; i < _b_values.size(); ++i )
  {
    int unique = 1;
    for(unsigned int j = 0; j < Bunique.size(); ++j )
    {
      if (_b_values[i] == Bunique[j])
      {
        unique = 0;
        break;
      }
    }
    if(unique == 1)
    {
      Bunique.push_back(_b_values[i]);
    }
    if(Bmax < _b_values[i])
    {
      Bmax = _b_values[i];
    }
  }

  tmp = cbrt(3*Bmax*1000000/(2*GAMMA*GAMMA*Gmax*Gmax));

  for(int i = 0; i < _b_values.size(); ++i )
  {
    l_pulseSeparation[i] = tmp;
  }

  for(unsigned int i = 0; i < Bunique.size(); ++i )
  {
    tmpG.push_back(std::sqrt(Bunique[i]/Bmax) * Gmax);
    // std::cout<< "\n tmpG:" << std::sqrt(Bunique[i]/Bmax) * Gmax;
  }

  for(unsigned int i = 0; i < Bunique.size(); ++i )
  {
    for(int j=0; j < _b_values.size(); j++)
    {
      if(_b_values[j] == Bunique[i])
      {
        l_gradientStrength[j] = tmpG[i];
      }
    }
  }
}

void Tractography::UnpackTensor(const ukfVectorType& b,       // b - bValues
                                const stdVec_t& u,            // u - directions
                                stdEigVec_t& s,   // s = signal values
                                stdEigVec_t& ret) // starting params [i][j] : i - signal number; j
                                                                        // - param
{

  // DEBUGGING
  // std::cout << "b's: ";
  // for (int i=0; i<b.size();++i) {
  //   std::cout << b[i] << ", ";
  // }

  assert(ret.size() == s.size() );

  // Build B matrix.
  const int signal_dim = _signal_data->GetSignalDimension();
/**
  * A special type for holding 6 elements of tensor for each signal
  *  Only used in tractography.cc
  */
  typedef Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, 6>  BMatrixType;
  BMatrixType B(signal_dim * 2, 6); //HACK: Eigen::Matrix<ukfPrecisionType, DYNAMIC, 6> ??

  for( int i = 0; i < signal_dim * 2; ++i )
    {
    const vec3_t & g = u[i];
    B(i, 0) = (-b[i]) * (g[0] * g[0]);
    B(i, 1) = (-b[i]) * (2.0 * g[0] * g[1]);
    B(i, 2) = (-b[i]) * (2.0 * g[0] * g[2]);
    B(i, 3) = (-b[i]) * (g[1] * g[1]);
    B(i, 4) = (-b[i]) * (2.0 * g[1] * g[2]);
    B(i, 5) = (-b[i]) * (g[2] * g[2]);
    }
  // Use QR decomposition to find the matrix representation of the tensor at the
  // seed point of the fiber. Raplacement of the gmm function gmm::least_squares_cg(..)
  Eigen::HouseholderQR<BMatrixType> qr(B);

  // Temporary variables.
  mat33_t D;

  if (this->debug) std::cout << "Estimating seed tensors:" << std::endl;

  // Unpack data
  for( stdEigVec_t::size_type i = 0; i < s.size(); ++i )
    {
    for( unsigned int j = 0; j < s[i].size(); ++j )
      {
      if( s[i][j] <= 0 )
        {
        s[i][j] = 10e-8;
        }
      s[i][j] = log(s[i][j]);
      }

    // The six tensor components.
    //TODO: this could be fixed size
    ukfVectorType d = qr.solve(s[i]);

    // symmetric diffusion tensor
    D(0, 0) = d[0];
    D(0, 1) = d[1];
    D(0, 2) = d[2];
    D(1, 0) = d[1];
    D(1, 1) = d[3];
    D(1, 2) = d[4];
    D(2, 0) = d[2];
    D(2, 1) = d[4];
    D(2, 2) = d[5];
    // Use singular value decomposition to extract the eigenvalues and the
    // rotation matrix (which gives us main direction of the tensor).
    // NOTE that svd can be used here only because D is a normal matrix

    // std::cout << "Tensor test: " << std::endl << D << std::endl;
   Eigen::JacobiSVD<ukfMatrixType> svd_decomp(D,Eigen::ComputeThinU );
   mat33_t Q = svd_decomp.matrixU();
   vec3_t sigma = svd_decomp.singularValues(); // diagonal() returns elements of a diag matrix as a vector.

    assert(sigma[0] >= sigma[1] && sigma[1] >= sigma[2]);
    if( Q.determinant()  < ukfZero )
      {
      Q = Q * (-ukfOne);
      }
    assert(Q.determinant() > ukfZero);

    // Extract the three Euler Angles from the rotation matrix.
    ukfPrecisionType phi, psi;
    const ukfPrecisionType theta = std::acos(Q(2, 2) );
    ukfPrecisionType epsilon = 1.0e-10;
    if( fabs(theta) > epsilon )
      {
      phi = atan2(Q(1, 2), Q(0, 2) );
      psi = atan2(Q(2, 1), -Q(2, 0) );
      }
    else
      {
      phi = atan2(-Q(0, 1), Q(1, 1) );
      psi = ukfZero;
      }

    ret[i].resize(9);
    ret[i][0] = Q(0, 0);
    ret[i][1] = Q(1, 0);
    ret[i][2] = Q(2, 0);
    ret[i][3] = theta;
    ret[i][4] = phi;
    ret[i][5] = psi;
    sigma = sigma * GLOBAL_TENSOR_PACK_VALUE; // NOTICE this scaling of eigenvalues. The values are scaled back in diffusion_euler()
    ret[i][6] = sigma[0];
    ret[i][7] = sigma[1];
    ret[i][8] = sigma[2];
    }
}

void Tractography::Follow3T(const int thread_id,
                            const size_t seed_index,
                            const SeedPointInfo& fiberStartSeed,
                            UKFFiber& fiber,
                            bool is_branching,
                            std::vector<SeedPointInfo>& branching_seeds,
                            std::vector<BranchingSeedAffiliation>& branching_seed_affiliation)
{
  int fiber_size = 100;
  int fiber_length = 0;
  assert(_model->signal_dim() == _signal_data->GetSignalDimension() * 2);

  // Unpack the fiberStartSeed information.
  vec3_t              x = fiberStartSeed.point;
  State              state = fiberStartSeed.state;
  ukfMatrixType    p(fiberStartSeed.covariance);
  ukfPrecisionType             fa = fiberStartSeed.fa;
  ukfPrecisionType             fa2 = fiberStartSeed.fa2;
  ukfPrecisionType             dNormMSE = 0; // no error at the fiberStartSeed
  ukfPrecisionType             trace = fiberStartSeed.trace;
  ukfPrecisionType             trace2 = fiberStartSeed.trace2;

  //  Reserving fiber array memory so as to avoid resizing at every step
  FiberReserve(fiber, fiber_size);

  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

  vec3_t m1 = fiberStartSeed.start_dir;
  vec3_t l1, m2, l2, m3, l3;

  // Tract the fiber.
  ukfMatrixType signal_tmp(_model->signal_dim(), 1);
  ukfMatrixType state_tmp(_model->state_dim(), 1);

  int stepnr = 0;
  while( true )
    {
    ++stepnr;

    Step3T(thread_id, x, m1, l1, m2, l2, m3, l3, fa, fa2, state, p, dNormMSE, trace, trace2);

    // Check if we should abort following this fiber. We abort if we reach the
    // CSF, if FA or GA get too small, if the curvature get's too high or if
    // the fiber gets too long.
    const bool is_brain = _signal_data->ScalarMaskValue(x) > 0; //_signal_data->Interp3ScalarMask(x) > 0.1;

    state_tmp.col(0) = state;
    _model->H(state_tmp, signal_tmp);

    const ukfPrecisionType mean_signal = s2adc(signal_tmp);
    const bool in_csf = (mean_signal < _mean_signal_min) ||
                        (fa < _fa_min);

    bool is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain || in_csf
        || stepnr > _max_length  // Stop if the fiber is too long
        || is_curving )
      {

      break;

      }

    if(fiber_length>=fiber_size)
      {
        // If fibersize is more than initally allocated size resizing further
        fiber_size += 100;
        FiberReserve(fiber, fiber_size);
      }

    if((stepnr+1)%_steps_per_record == 0)
      {
        fiber_length++;
        Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);
      }

    // Record branch if necessary.
    if( is_branching )
      {
      const ukfPrecisionType fa_1 = l2fa(l1[0], l1[1], l1[2]);
      const bool is_one = ( l1[0] > l1[1] ) && (l1[0] > l1[2]) && ( fa_1 > _fa_min );
      if( is_one )
        {

        bool add_m2 = false;
        bool add_m3 = false;

        const ukfPrecisionType fa_2 = l2fa(l2[0], l2[1], l2[2]);
        const ukfPrecisionType fa_3 = l2fa(l3[0], l3[1], l3[2]);

        const bool is_two = ( l2[0] > l2[1] ) && ( l2[0] > l2[2]) && ( fa_2 > _fa_min );
        const bool is_three = ( l3[0] > l3[1] ) && ( l3[0] > l3[2] ) && ( fa_3 > _fa_min );

        ukfPrecisionType dotval = m1.dot(m2);
        const bool is_branch1 = dotval < _cos_theta_min && dotval > _cos_theta_max;
        dotval = m1.dot(m3);
        const bool is_branch2 = dotval < _cos_theta_min && dotval > _cos_theta_max;
        dotval = m2.dot(m3);
        const bool is_branch3 = dotval < _cos_theta_min;

        int state_dim = _model->state_dim();
        // If there is a branch between m1 and m2.
        if( is_two && is_branch1 )
          {
          // If there is a branch between m1 and m3 we have to check if we
          // branch twice or only once.
          if( is_three && is_branch2 )
            {
            // If angle between m2 and m3 is big enough we have 2 branches.
            if( is_branch3 )
              {
              add_m2 = true;
              add_m3 = true;
              }
            else
              {
              // Otherwise we only follow m2 or m3, and we follow the one
              // tensor where the FA is bigger.
              if( fa_2 > fa_3 )
                {
                add_m2 = true;
                }
              else
                {
                add_m3 = true;
                }
              }
            }
          else
            {
            // If it's not possible for m3 to branch we are sure that m2
            // branches.
            add_m2 = true;
            }
          }
        else if( is_three && is_branch2 )
          {
          // If m2 is not branching we only check if m3 can branch.
          add_m3 = true;
          }

        // If we have to tensors and the angle between them is large enough we
        // create a new seed for the branch. Since this branch is following the
        // second tensor we swap the state and covariance.
        if( add_m2 )
          {
          branching_seeds.push_back(SeedPointInfo() );
          SeedPointInfo& local_seed = branching_seeds[branching_seeds.size() - 1];
          branching_seed_affiliation.push_back(BranchingSeedAffiliation() );
          BranchingSeedAffiliation& affiliation = branching_seed_affiliation[branching_seed_affiliation.size() - 1];

          affiliation.fiber_index_ = seed_index;
          affiliation.position_on_fiber_ = stepnr;

          local_seed.state.resize(state_dim);
          local_seed.state = state;

          local_seed.covariance.resize(state_dim, state_dim);
          local_seed.covariance = p;

          SwapState3T(local_seed.state, local_seed.covariance, 2);
          local_seed.point = x;
          local_seed.start_dir = m2;
          local_seed.fa = fa_2;
          }
        // Same for the third tensor.
        if( add_m3 )
          {
          branching_seeds.push_back(SeedPointInfo() );
          SeedPointInfo& local_seed = branching_seeds[branching_seeds.size() - 1];
          branching_seed_affiliation.push_back(BranchingSeedAffiliation() );
          BranchingSeedAffiliation& affiliation = branching_seed_affiliation[branching_seed_affiliation.size() - 1];

          affiliation.fiber_index_ = seed_index;
          affiliation.position_on_fiber_ = stepnr;

          local_seed.state.resize(state_dim);
          local_seed.state = state;
          local_seed.covariance.resize(state_dim, state_dim);
          local_seed.covariance = p;
          SwapState3T(local_seed.state, local_seed.covariance, 3);
          local_seed.point = x;
          local_seed.start_dir = m3;
          local_seed.fa = fa_3;
          }
        }
      }
    }
  FiberReserve(fiber, fiber_length);
}

void Tractography::Follow2T(const int thread_id,
                            const size_t seed_index,
                            const SeedPointInfo& fiberStartSeed,
                            UKFFiber& fiber,
                            bool is_branching,
                            std::vector<SeedPointInfo>& branching_seeds,
                            std::vector<BranchingSeedAffiliation>& branching_seed_affiliation)
{
  int fiber_size = 100;
  int fiber_length = 0;
  // Unpack the fiberStartSeed information.
  vec3_t x = fiberStartSeed.point;   // NOTICE that the x here is in ijk coordinate system
  State state = fiberStartSeed.state;

  ukfMatrixType p(fiberStartSeed.covariance);

  ukfPrecisionType fa = fiberStartSeed.fa;
  ukfPrecisionType fa2 = fiberStartSeed.fa2;
  ukfPrecisionType dNormMSE = 0; // no error at the fiberStartSeed
  ukfPrecisionType trace = fiberStartSeed.trace;
  ukfPrecisionType trace2 = fiberStartSeed.trace2;

  //  Reserving fiber array memory so as to avoid resizing at every step
  FiberReserve(fiber, fiber_size);

  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2); // writes state to fiber.state

  vec3_t m1, l1, m2, l2;
  m1 = fiberStartSeed.start_dir;

  // Track the fiber.
  ukfMatrixType signal_tmp(_model->signal_dim(), 1);
  ukfMatrixType state_tmp(_model->state_dim(), 1);

  int stepnr = 0;

  // useful for debuggingo
//   std::ofstream stateFile;
//   stateFile.open("states.txt", std::ios::app);

  while( true )
    {
    ++stepnr;

    Step2T(thread_id, x, m1, l1, m2, l2, fa, fa2, state, p, dNormMSE, trace, trace2);

    // Check if we should abort following this fiber. We abort if we reach the
    // CSF, if FA or GA get too small, if the curvature get's too high or if
    // the fiber gets too long.
    const bool is_brain = _signal_data->ScalarMaskValue(x) > 0; // _signal_data->Interp3ScalarMask(x) > 0.1; // is this 0.1 correct? yes

    // after here state does not change until next step.

    state_tmp.col(0) = state;

    _model->H(state_tmp, signal_tmp); // signal_tmp is written, but only used to calculate mean signal

    const ukfPrecisionType mean_signal = s2adc(signal_tmp);
    const bool in_csf = (_noddi) ? ( mean_signal < _mean_signal_min ) :
                                   ( mean_signal < _mean_signal_min || fa < _fa_min);

    const bool is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain
        || in_csf
        || stepnr > _max_length  // Stop when the fiber is too long
        || is_curving)
      {
      break;
      }

    if (_noddi)
      if (state[4] < 0.6 || state[9] < 0.6) // kappa1 and kappa2 break conditions
        break;

    if(fiber_length>=fiber_size)
      {
        // If fibersize is more than initally allocated size resizing further
        fiber_size += 100;
        FiberReserve(fiber, fiber_size);
      }

    if((stepnr+1)%_steps_per_record == 0)
      {
        fiber_length++;
        if(_noddi)
          Record(x, state[3], state[8], state, p, fiber, dNormMSE, state[4], state[9]);
        else
          Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);
      }

    // Record branch if necessary.
    if( is_branching )
      {
      const ukfPrecisionType fa_1 = l2fa(l1[0], l1[1], l1[2]);
      const ukfPrecisionType fa_2 = l2fa(l2[0], l2[1], l2[2]);
      const bool is_two = ( l1[0] > l1[1] ) && ( l1[0] > l1[2] ) && ( l2[0] > l2[1] ) && (l2[0] > l2[2] )
        && ( fa_1 > _fa_min )&& ( fa_2 > _fa_min );
      ukfPrecisionType theta = m1.dot(m2);
      bool   is_branch = theta<_cos_theta_min && theta> _cos_theta_max;

      // If we have two tensors and the angle between them is large enough we
      // create a new seed for the branch. Since this branch is following the
      // second tensor we swap the state and covariance.
      if( is_two && is_branch )
        {
        branching_seeds.push_back(SeedPointInfo() );
        SeedPointInfo& local_seed = branching_seeds[branching_seeds.size() - 1];
        branching_seed_affiliation.push_back(BranchingSeedAffiliation() );
        BranchingSeedAffiliation& affiliation = branching_seed_affiliation[branching_seed_affiliation.size() - 1];

        affiliation.fiber_index_ = seed_index;
        affiliation.position_on_fiber_ = stepnr;

        int state_dim = _model->state_dim();
        local_seed.state.resize(state_dim);
        local_seed.state = state;
        local_seed.covariance.resize(state_dim, state_dim);
        local_seed.covariance = p;
        SwapState2T(local_seed.state, local_seed.covariance);
        local_seed.point = x;
        local_seed.start_dir = m2;
        local_seed.fa = fa_2;
        }
      }
    }
    FiberReserve(fiber, fiber_length);
//   stateFile.close();
}

// Also read the comments to Follow2T above, it's documented better than this
// function here.
void Tractography::Follow1T(const int thread_id,
                            const SeedPointInfo& fiberStartSeed,
                            UKFFiber& fiber)
{
  int fiber_size = 100;
  int fiber_length = 0;
  assert(_model->signal_dim() == _signal_data->GetSignalDimension() * 2);

  vec3_t x = fiberStartSeed.point;
  State state = fiberStartSeed.state;

  // DEBUG
//   std::cout << "fiberStartSeed state:\n";
//   for (int i=0;i<state.size();++i) {
//     std::cout << state[i] << " ";
//   }
//   std::cout << std::endl;

  ukfMatrixType p(fiberStartSeed.covariance);

  ukfPrecisionType fa = fiberStartSeed.fa;
  ukfPrecisionType fa2 = fiberStartSeed.fa2; // just needed for record
  ukfPrecisionType trace = fiberStartSeed.trace;
  ukfPrecisionType trace2 = fiberStartSeed.trace2;

  ukfPrecisionType dNormMSE = 0; // no error at the fiberStartSeed

  //  Reserving fiber array memory so as to avoid resizing at every step
  FiberReserve(fiber, fiber_size);

  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

  // Tract the fiber.
  ukfMatrixType signal_tmp(_model->signal_dim(), 1);
  ukfMatrixType state_tmp(_model->state_dim(), 1);

  int stepnr = 0;
  while( true )
    {
    ++stepnr;

    Step1T(thread_id, x, fa, state, p, dNormMSE, trace);

    // Terminate if off brain or in CSF.
    const bool is_brain = _signal_data->ScalarMaskValue(x) > 0; //_signal_data->Interp3ScalarMask(x) > 0.1; // x is the seed point
    state_tmp.col(0) = state;

    _model->H(state_tmp, signal_tmp);

    // Check mean_signal threshold
    const ukfPrecisionType mean_signal = s2adc(signal_tmp);
    bool in_csf;
    if(_noddi)
      in_csf = mean_signal < _mean_signal_min;
    else
      in_csf = mean_signal < _mean_signal_min || fa < _fa_min;

    bool is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain
        || in_csf
        || stepnr > _max_length  // Stop when fiber is too long
        || is_curving )
      {
      break;
      }
    if (_noddi)
      if( state[4]<1.2) // checking kappa
        break;

    if(fiber_length>=fiber_size)
      {
        // If fibersize is more than initally allocated size resizing further
        fiber_size += 100;
        FiberReserve(fiber, fiber_size);
      }

    if((stepnr+1)%_steps_per_record == 0)
      {
        fiber_length++;
        if(_noddi)
          Record(x, state[3], fa2, state, p, fiber, dNormMSE, state[4], trace2);
        else
          Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);
      }


    }
    FiberReserve(fiber, fiber_length);
}

void Tractography::Step3T(const int thread_id,
                          vec3_t& x,
                          vec3_t& m1,
                          vec3_t& l1,
                          vec3_t& m2,
                          vec3_t& l2,
                          vec3_t& m3,
                          vec3_t& l3,
                          ukfPrecisionType& fa,
                          ukfPrecisionType& fa2,
                          State& state,
                          ukfMatrixType& covariance,
                          ukfPrecisionType& dNormMSE,
                          ukfPrecisionType& trace,
                          ukfPrecisionType& trace2
                          )
{

  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
         static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );
  State state_new(_model->state_dim() );

  ukfMatrixType covariance_new(_model->state_dim(), _model->state_dim() );

  // Use the Unscented Kalman Filter to get the next estimate.
  ukfVectorType signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);
  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  vec3_t old_dir = m1;

  _model->State2Tensor3T(state, old_dir, m1, l1, m2, l2, m3, l3);
  trace = l1[0] + l1[1] + l1[2];
  trace2 = l2[0] + l2[1] + l2[2];

  ukfPrecisionType dot1 = m1.dot(old_dir);
  ukfPrecisionType dot2 = m2.dot(old_dir);
  ukfPrecisionType dot3 = m3.dot(old_dir);
  if( dot1 < dot2 && dot3 < dot2 )
    {
    // Switch dirs and lambdas.
    vec3_t tmp = m1;
    m1 = m2;
    m2 = tmp;
    tmp = l1;
    l1 = l2;
    l2 = tmp;

    // Swap state.

    SwapState3T(state, covariance, 2);

    }
  else if( dot1 < dot3 )
    {
    // Switch dirs and lambdas.
    vec3_t tmp = m1;
    m1 = m3;
    m3 = tmp;
    tmp = l1;
    l1 = l3;
    l3 = tmp;

    // Swap state.
    SwapState3T(state, covariance, 3);
    }

  // Update FA. If the first lamba is not the largest anymore the FA is set to
  // 0, and the 0 FA value will lead to abortion in the tractography loop.
  if( l1[0] < l1[1] || l1[0] < l1[2] )
    {
    fa = ukfZero;
    }
  else
    {
    fa = l2fa(l1[0], l1[1], l1[2]);
    fa2 = l2fa(l2[0], l2[1], l2[2]);
    }

  const vec3_t & voxel = _signal_data->voxel();

  // CB: Bug corrected, dir[i] should be divided by voxel[i]
  vec3_t dx;
  dx << m1[2] / voxel[0],
    m1[1] / voxel[1],
    m1[0] / voxel[2];
  x = x + dx * _stepLength;
}

void Tractography::Step2T(const int thread_id,
                          vec3_t& x,
                          vec3_t& m1,
                          vec3_t& l1,
                          vec3_t& m2,
                          vec3_t& l2,
                          ukfPrecisionType& fa,
                          ukfPrecisionType& fa2,
                          State& state,
                          ukfMatrixType& covariance,
                          ukfPrecisionType& dNormMSE,
                          ukfPrecisionType& trace,
                          ukfPrecisionType& trace2
                          )
{
  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
    static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );

  State              state_new(_model->state_dim() );
  ukfMatrixType covariance_new(_model->state_dim(), _model->state_dim() );
  covariance_new.setConstant(ukfZero);

  // Use the Unscented Kalman Filter to get the next estimate.
  ukfVectorType signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);

  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  const vec3_t old_dir = m1;   // Direction in last step
  ukfPrecisionType fa_tensor_1 = ukfZero;
  ukfPrecisionType fa_tensor_2 = ukfZero;
  if(_noddi)
    {
    initNormalized(m1, state[0], state[1], state[2]);
    initNormalized(m2, state[5], state[6], state[7]);
    if( m1[0] * old_dir[0] + m1[1] * old_dir[1] + m1[2] * old_dir[2] < 0 )
      {
      m1 = -m1;
      }
    if( m2[0] * old_dir[0] + m2[1] * old_dir[1] + m2[2] * old_dir[2] < 0 )
      {
      m2 = -m2;
      }
    }
  else
    {
    _model->State2Tensor2T(state, old_dir, m1, l1, m2, l2);   // The returned m1 and m2 are unit vector here
    trace = l1[0] + l1[1] + l1[2];
    trace2 = l2[0] + l2[1] + l2[2];
    fa_tensor_1 = l2fa(l1[0], l1[1], l1[2]);
    fa_tensor_2 = l2fa(l2[0], l2[1], l2[2]);
    }
  const ukfPrecisionType tensor_angle = RadToDeg( std::acos(m1.dot(m2) ) );

  if( m1.dot(old_dir) < m2.dot(old_dir) )
    {
    // Switch dirs and lambdas.
    vec3_t tmp = m1;
    m1 = m2;
    m2 = tmp;
    tmp = l1;
    l1 = l2;
    l2 = tmp;
    // Need to swap scalar measures too.
    ukfPrecisionType tmpScalar = fa_tensor_1;
    fa_tensor_1 = fa_tensor_2;
    fa_tensor_2=tmpScalar;
    tmpScalar= trace;
    trace= trace2;
    trace2=tmpScalar;
    ukfMatrixType old = covariance;
    SwapState2T(state, covariance);   // Swap the two tensors
    }

  if(tensor_angle <= 20)
    {
    if(_noddi)
      {
      vec3_t tmp = m1;
      m1 = m2;
      m2 = tmp;
      ukfMatrixType old = covariance;

      SwapState2T(state, covariance);   // Swap the two tensors
      }
    else if (std::min(fa_tensor_1, fa_tensor_2) <= 0.2 )
      {
      if( fa_tensor_1 > 0.2 )
        {
        // do nothing
        // i.e. keep m1 as principal direction
        }
      else
        {
        // switch directions, note: FA will be re-calculated after
        vec3_t tmp = m1;
        m1 = m2;
        m2 = tmp;
        tmp = l1;
        l1 = l2;
        l2 = tmp;

        ukfMatrixType old = covariance;

        SwapState2T(state, covariance);   // Swap the two tensors
        }
      }
    }
  // Update FA. If the first lamba is not the largest anymore the FA is set to
  // 0 what will lead to abortion in the tractography loop.
  if( l1[0] < l1[1] || l1[0] < l1[2] )
    {
    fa = ukfZero;
    fa2 = ukfZero;
    }
  else
    {
    fa = l2fa(l1[0], l1[1], l1[2]);
    fa2 = l2fa(l2[0], l2[1], l2[2]);
    }
  vec3_t dx;
    {
    const vec3_t dir = m1;     // The dir is a unit vector in ijk coordinate system indicating the direction of step
    const vec3_t voxel = _signal_data->voxel();
    dx << dir[2] / voxel[0], // By dividing by the voxel size, it's guaranteed that the step
       // represented by dx is 1mm in RAS coordinate system, no matter whether
       // the voxel is isotropic or not
       dir[1] / voxel[1], // The value is scaled back during the ijk->RAS transformation when
       // outputted
       dir[0] / voxel[2];

    x = x + dx * _stepLength; // The x here is in ijk coordinate system.
    }
  // NOTICE that the coordinate order of x is in reverse order with respect to the axis order in the original signal
  // file.
  // This coordinate order is filpped back during output
  // The step length is in World space
  // throw;
}

void Tractography::Step1T(const int thread_id,
                          vec3_t& x,
                          ukfPrecisionType& fa,
                          State& state,
                          ukfMatrixType& covariance,
                          ukfPrecisionType& dNormMSE,
                          ukfPrecisionType& trace
                          )
{

  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
         static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );
  State              state_new(_model->state_dim() );
  ukfMatrixType covariance_new(_model->state_dim(), _model->state_dim() );

  ukfVectorType signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);

  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  vec3_t dir;
  if (_noddi)
  {
    dir << state[0], state[1], state[2];
  }
  else
  {
    vec3_t l;
    _model->State2Tensor1T(state, dir, l);

    trace = l[0] + l[1] + l[2];

    // Update FA. If the first lamba is not the largest anymore the FA is set to
    // 0 what will lead to abortion in the tractography loop.
    if( l[0] < l[1] || l[0] < l[2] )
      {
      fa = ukfZero;
      }
    else
      {
      fa = l2fa(l[0], l[1], l[2]);
      }
  }
  vec3_t voxel = _signal_data->voxel();

  vec3_t dx;
  dx << dir[2] / voxel[0],
    dir[1] / voxel[1],
    dir[0] / voxel[2];
  x = x + dx * _stepLength;

}

void Tractography::SwapState3T(State& state,
                               ukfMatrixType& covariance,
                               int i)
{

  // This function is only for 3T.
  assert(i == 2 || i == 3);


  int                state_dim = _model->state_dim();
  ukfMatrixType tmp(state_dim, state_dim);
  state_dim /= 3;
  assert(state_dim == 5 || state_dim == 6);
  --i;
  int j = i == 1 ? 2 : 1;
  i *= state_dim;
  j *= state_dim;

  tmp.setConstant(ukfZero);
  tmp = covariance;
  covariance.block( i, i,state_dim, state_dim) = tmp.block( 0, 0,state_dim, state_dim);
  covariance.block( 0, 0,state_dim, state_dim) = tmp.block( i, i,state_dim, state_dim);
  covariance.block( 0, i,state_dim, state_dim) = tmp.block( i, 0,state_dim, state_dim);
  covariance.block( i, 0,state_dim, state_dim) = tmp.block( 0, i,state_dim, state_dim);

  covariance.block( j, i,state_dim, state_dim) = tmp.block( j, 0,state_dim, state_dim);
  covariance.block( j, 0,state_dim, state_dim) = tmp.block( j, i,state_dim, state_dim);
  covariance.block( i, j,state_dim, state_dim) = tmp.block( 0, j,state_dim, state_dim);
  covariance.block( 0, j,state_dim, state_dim) = tmp.block( i, j,state_dim, state_dim);

  // Swap the state.
  const ukfVectorType tmp_vec = state;
  state.segment(i,state_dim) = tmp_vec.segment(0,state_dim);
  state.segment(0,state_dim) = tmp_vec.segment(i,state_dim);
}

void Tractography::SwapState2T( State& state,
                                ukfMatrixType& covariance)
{
  // This function is only for 2T.
  int state_dim = _model->state_dim();

  ukfMatrixType tmp(state_dim, state_dim);
  bool               bUnevenState = false;

  if( state_dim % 2 != 0 )
    {
    bUnevenState = true;                       // there is a weight term in the end of the state
    }
  state_dim = state_dim >> 1; // for uneven state (fw) rounds down, thats good

  tmp.setConstant(ukfZero);
  tmp = covariance;

  covariance.block( state_dim, state_dim,state_dim, state_dim) = tmp.block( 0, 0,state_dim, state_dim);
  covariance.block( 0, 0,state_dim, state_dim) = tmp.block( state_dim, state_dim,state_dim, state_dim);
  covariance.block( 0, state_dim,state_dim, state_dim) = tmp.block( state_dim, 0,state_dim, state_dim);
  covariance.block( state_dim, 0,state_dim, state_dim) = tmp.block( 0, state_dim,state_dim, state_dim);

  if( bUnevenState )   // change covariances of weights and state so they match the state again
    {
    covariance.block( state_dim * 2, state_dim,1, state_dim) = tmp.block( state_dim * 2, 0,1, state_dim);
    covariance.block( state_dim * 2, 0,1, state_dim) = tmp.block( state_dim * 2, state_dim,1, state_dim);

    covariance.block( state_dim, state_dim * 2,state_dim, 1) = tmp.block( 0, state_dim * 2,state_dim, 1);
    covariance.block( 0, state_dim * 2,state_dim, 1) = tmp.block( state_dim, state_dim * 2,state_dim, 1);
    }

  // Swap the state.
  const ukfVectorType tmp_vec = state;
  state.segment(state_dim,state_dim) = tmp_vec.segment(0,state_dim);
  state.segment(0,state_dim) = tmp_vec.segment(state_dim, state_dim);
}

void Tractography::Record(const vec3_t& x, const ukfPrecisionType fa, const ukfPrecisionType fa2, const State& state,
                          const ukfMatrixType p,
                          UKFFiber& fiber, const ukfPrecisionType dNormMSE, const ukfPrecisionType trace, const ukfPrecisionType trace2)
{
  // if Noddi model is used Kappa is stored in trace, Vic in fa and Viso in freewater
  assert(_model->state_dim() == static_cast<int>(state.size() ) );
  assert(p.rows() == static_cast<unsigned int>(state.size() ) &&
         p.cols() == static_cast<unsigned int>(state.size() ) );

  // std::cout << "x: " << x[0] << " " << x[1] << " " << x[2] << std::endl;
  fiber.position.push_back(x);
  fiber.norm.push_back(p.norm());

  if( _record_nmse )
    {
    fiber.normMSE.push_back(dNormMSE);
    }

  if( _record_trace || _record_kappa)
    {
    fiber.trace.push_back(2*(atan(1/trace)/3.14));
    if( _num_tensors >= 2 )
      {
      fiber.trace2.push_back(2*(atan(1/trace2)/3.14));
      }
    }

  if( _record_fa || _record_Vic)
    {
    fiber.fa.push_back(fa);
    if( _num_tensors >= 2 )
      {
      fiber.fa2.push_back(fa2);
      }
    }

  if ( _record_Viso )
    {
    ukfPrecisionType viso = state[_nPosFreeWater];
    if( viso < 0 )
      {
      if( viso >= -1.0e-4 ) // for small errors just round it to 0
        {
        viso = 0;
        }
      else   // for too big errors exit with exception.
        {
        std::cout << "Error: program produced negative free water.\n";
        throw;
        }
      }
    fiber.free_water.push_back(viso);
    }

  if( _record_free_water )
    {
    ukfPrecisionType fw = 1 - state[_nPosFreeWater];
    // sometimes QP produces slightly negative results due to numerical errors in Quadratic Programming, the weight is
    // clipped in F() and H() but its still possible that
    // a slightly negative weight gets here, because the filter ends with a constrain step.
    if( fw < 0 )
      {
      if( fw >= -1.0e-4 ) // for small errors just round it to 0
        {
        fw = 0;
        }
      else   // for too big errors exit with exception.
        {
        std::cout << "Error: program produced negative free water.\n";
        throw;
        }
      }
    fiber.free_water.push_back(fw);
    }

  // Record the state
  if( state.size() == 5 || state.size() == 10 || state.size() == 15 )   // i.e. simple model
    { // Normalize direction before storing it;
    State store_state(state);
    vec3_t dir;
    initNormalized(dir, store_state[0], store_state[1], store_state[2]);
    store_state[0] = dir[0];
    store_state[1] = dir[1];
    store_state[2] = dir[2];

    if( state.size() == 10 )
      {
      initNormalized(dir,store_state[5], store_state[6], store_state[7]);
      store_state[5] = dir[0];
      store_state[6] = dir[1];
      store_state[7] = dir[2];
      }
    if( state.size() == 15 )
      {
      initNormalized(dir,store_state[10], store_state[11], store_state[12]);
      store_state[10] = dir[0];
      store_state[11] = dir[1];
      store_state[12] = dir[2];
      }
    fiber.state.push_back(store_state);

    }
  else
    {
    // Normal state
    fiber.state.push_back(state);
    }

  if( _record_cov )
    {
    fiber.covariance.push_back(p);
    }
}


void Tractography::FiberReserve(UKFFiber& fiber, int fiber_size)
{
  // Reserving space for fiber
  fiber.position.reserve(fiber_size);
  fiber.norm.reserve(fiber_size);
  fiber.state.reserve(fiber_size);
  if( _record_nmse )
    {
    fiber.normMSE.reserve(fiber_size);
    }

  if( _record_trace || _record_kappa)
    {
    fiber.trace.reserve(fiber_size);
    if( _num_tensors >= 2 )
      {
      fiber.trace2.reserve(fiber_size);
      }
    }

  if( _record_fa || _record_Vic)
    {
    fiber.fa.reserve(fiber_size);
    if( _num_tensors >= 2 )
      {
      fiber.fa2.reserve(fiber_size);
      }
    }
  if( _record_free_water || _record_Viso)
    {
    fiber.free_water.reserve(fiber_size);
    }
  if( _record_cov)
    {
      fiber.covariance.reserve(fiber_size);
    }
}
