/**
 * \file tractography.cc
 * \brief implementation of tractography.h
*/

#include "tractography.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
// #include "timer.h"
#include "filter_model.h"
#include "ISignalData.h"
#include "NrrdData.h"
#include "utilities.h"
#include "vtk_writer.h"
#include "thread.h"
// #include "UKFMultiThreader.h"
#include "itkMultiThreader.h"
#include "math_utilities.h"
#include "filter_Simple2BiExp_FW.h"

#include <fstream>
#include <iostream>

// TODO implement this switch
#include "config.h"

Tractography::Tractography(FilterModel *model, model_type filter_model_type,

                           const std::string& output_file, const std::string & output_file_with_second_tensor,
                           const bool record_fa, const bool record_nmse, const bool record_trace,
                           const bool record_state,
                           const bool record_cov, const bool record_free_water,  const bool record_tensors,
                           const bool record_Vic, const bool record_kappa,  const bool record_Viso,
                           const bool transform_position, const bool store_glyphs, const bool branchesOnly,

                           const ukfPrecisionType fa_min, const ukfPrecisionType ga_min, const ukfPrecisionType seedFALimit,
                           const int num_tensors, const int seeds_per_voxel,
                           const ukfPrecisionType minBranchingAngle, const ukfPrecisionType maxBranchingAngle,
                           const bool is_full_model, const bool free_water, const bool noddi,
                           const bool diffusionPropagator, const ukfPrecisionType rtop_min, const bool recordRTOP, 
                           const ukfPrecisionType maxNMSE, const ukfPrecisionType maxUKFIterations,
   
                           const ukfPrecisionType stepLength, const ukfPrecisionType recordLength,
                           const ukfPrecisionType maxHalfFiberLength,
                           const std::vector<int>& labels,

                           ukfPrecisionType p0, ukfPrecisionType sigma_signal, ukfPrecisionType sigma_mask,
                           ukfPrecisionType min_radius,
                           ukfPrecisionType /* UNUSED full_brain_ga_min */,

                           const int num_threads
                           ) :
  _ukf(0, NULL), _model(model), _filter_model_type(filter_model_type),

  _output_file(output_file), _output_file_with_second_tensor(output_file_with_second_tensor),
  _record_fa(record_fa), _record_nmse(record_nmse), _record_trace(record_trace), _record_state(record_state),
  _record_cov(record_cov), _record_free_water(record_free_water),
  _record_Vic(record_Vic),
  _record_kappa(record_kappa), _record_Viso(record_Viso),
  _record_tensors(record_tensors),
  _transform_position(transform_position), _store_glyphs(store_glyphs), _branches_only(branchesOnly),

  _p0(p0), _sigma_signal(sigma_signal), _sigma_mask(sigma_mask), _min_radius(min_radius),
  //UNUSED _full_brain_ga_min(full_brain_ga_min),
  _max_length(static_cast<int>(std::ceil(maxHalfFiberLength / stepLength) ) ),
  _full_brain(false),
  _noddi(noddi),
  _fa_min(fa_min), _ga_min(ga_min), _seedFALimit(seedFALimit),
  _num_tensors(num_tensors), _seeds_per_voxel(seeds_per_voxel),
  _cos_theta_min(minBranchingAngle), _cos_theta_max(maxBranchingAngle),
  _is_full_model(is_full_model),
  _diffusionPropagator(diffusionPropagator), _rtop_min(rtop_min), _recordRTOP(recordRTOP), 
  _maxNMSE(maxNMSE), _maxUKFIterations(maxUKFIterations),
  _free_water(free_water),
  _stepLength(stepLength),
  _steps_per_record(recordLength/stepLength),
  _labels(labels),
  _writeBinary(true),
  _writeCompressed(true),
  _num_threads(num_threads)
{
  if( _cos_theta_max != ukfZero && _cos_theta_max <= _cos_theta_min )
    {
    std::cout << "Maximum branching angle must be greater than " << minBranchingAngle << " degrees." << std::endl;
    exit(1);
    }

  if( _num_tensors < 1 || _num_tensors > 3 )
    {
    std::cout << "Only one, two or three tensors are supported." << std::endl;
    exit(1);
    }

  _cos_theta_max = std::cos( DegToRad( _cos_theta_max ) );
  _cos_theta_min = std::cos( DegToRad( _cos_theta_min ) );

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
    if( !is_full_model && free_water ) // simple model with free water
      {
      _nPosFreeWater = 5;
      }
    else if( is_full_model && free_water ) // full model with free water
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
    if (_diffusionPropagator) {
      _nPosFreeWater = 14;
     }
    else if( !is_full_model && free_water ) // simple model with free water
      {
      _nPosFreeWater = 10;
      }
    else if( is_full_model && free_water ) // full model with free water
      {
      _nPosFreeWater = 12;
      }
    else if( _noddi ) // Viso is recorded
      {
      _nPosFreeWater = 10;
      }
    }
}

Tractography::~Tractography()
{
  if( _signal_data )
    {
    delete _signal_data;
    }
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

  _model->set_signal_data(_signal_data);

  _model->set_signal_dim(_signal_data->GetSignalDimension() * 2);

  return false;
}

void Tractography::Init(std::vector<SeedPointInfo>& seed_infos)
{

  assert(_signal_data);
  int signal_dim = _signal_data->GetSignalDimension();

  stdVec_t seeds;
  assert(_labels.size() > 0);
  if( !_full_brain )
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
          if( _signal_data->Interp3ScalarMask(pos) > 0.1 )
            {
            seeds.push_back(pos);
            }
          }
        }
      }
    }

  assert(seeds.size() > 0);

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
  int num_ga_too_low = 0;

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
        if( _full_brain && s2ga(signal_tmp) < _seedFALimit )
          {
          keep = false;
          ++num_ga_too_low;
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
  std::cout<<"Processing "<<starting_points.size()<<" starting points"<<std::endl;
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

    if( fa <= _seedFALimit )
      {
      ++fa_too_low;
      continue;
      }
    // In the case of the diffusion propagator model, if lambda1 < lambda2, 
    // then we do not consider this seed
    
    if (_diffusionPropagator && param[6] < param[7]) {
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
    } // If diffusion propagator
    else if (_diffusionPropagator) {
    
      //
      // STEP 1: Initialise the state based on the single estimated tensor
      //
      tmp_info_state.resize(15);
      tmp_info_inv_state.resize(15);
      
      // Diffusion directions, m1 = m2
      tmp_info_state[0] = tmp_info_state[7] = info.start_dir[0];
      tmp_info_state[1] = tmp_info_state[8] = info.start_dir[1];
      tmp_info_state[2] = tmp_info_state[9] = info.start_dir[2];
      
      // Fast diffusing component,  lambdas l11, l21 = l1 from the single tensor
      //                            lambdas l12, l21 = (l2 + l3) /2 from the single tensor (avg already calculated and stored in l2)
      tmp_info_state[3] = tmp_info_state[10] = param[6];
      tmp_info_state[4] = tmp_info_state[11] = param[7];
      
      // Slow diffusing component,  lambdas l13, l23 = 0.2 * l1
      //                            lambdas l14, l24 = 0.2 * (l2 + l3) /2
      // tmp_info_state[5] = tmp_info_state[12] = 0.2 * param[6];
      // tmp_info_state[6] = tmp_info_state[13] = 0.2 * param[7];
      tmp_info_state[5] = tmp_info_state[12] = 0.7 * param[6];
      tmp_info_state[6] = tmp_info_state[13] = 0.7 * param[7];
      
      // Free water volume fraction
      tmp_info_state[14] = 0.9; // 0.9 as initial value

      //
      // STEP 2.1: Loop the UKF at the same point in space. 
      // The UKF is an estimator, and we want to find the estimate with the smallest error through the iterations
      
      // Set the covariance value
      const int state_dim = tmp_info_state.size();
      info.covariance.resize(state_dim, state_dim);
      info_inv.covariance.resize(state_dim, state_dim);
      
      // make sure covariances are really empty
      info.covariance.setConstant(ukfZero);
      info_inv.covariance.setConstant(ukfZero);
      
      // fill the diagonal of the covariance matrix with _p0 (zeros elsewhere)
      for (int local_i = 0; local_i < state_dim; ++local_i) {
        info.covariance(local_i, local_i) = _p0;
        info_inv.covariance(local_i, local_i) = _p0;
      }
      
      // Input of the filter
      State state = ConvertVector<stdVecState, State>(tmp_info_state);
      ukfMatrixType p(info.covariance);
      
      ukfPrecisionType dNormMSE = 0.0;
      
      // Estimate the initial state
      
      //InitLoopUKF(state, p, signal_values[i], dNormMSE);
      NonLinearLeastSquareOptimization(state, signal_values[i], _model);

      
      
      // Output of the filter
      tmp_info_state = ConvertVector<State, stdVecState>(state);
     
      //
      // STEP 2.2: Filter out the seeds with a dNMSE > _maxNMSE (0.15 by default)
      if (dNormMSE > _maxNMSE) {
        continue;
      }
      
      //
      // STEP 3: Compute the return to origin probability (RTOP)
      // The RTOP is computed using the state informations
      ukfPrecisionType rtop = 0.0;
      ukfPrecisionType rtop1 = 0.0;
      ukfPrecisionType rtop2 = 0.0;
      ukfPrecisionType rtopSignal = 0.0;
      
      computeRTOPfromState(tmp_info_state, rtop, rtop1, rtop2); 
      computeRTOPfromSignal(rtopSignal, signal_values[i]);
      
      // These values are stored so that: rtop1 -> fa; rtop2 -> fa2; rtop -> trace; rtopSignal -> trace2
      
      info.fa = rtop1;
      info.fa2 = rtop2;
      info.trace = rtop;
      info.trace2 = rtopSignal;
      
      info_inv.fa = rtop1;
      info_inv.fa2 = rtop2;
      info_inv.trace = rtop;
      info_inv.trace2 = rtopSignal;
      
      //
      // STEP 4: Adjust the seed state info, based on the number of branches found at the seed point.
      //
      
      // Angle between m1 and m2 (m1 and m2 are unit vectors)
      
      
      const ukfPrecisionType min_angle = 30; // 25 or 30
      ukfPrecisionType scalar_product_m1_m2 = tmp_info_state[0] * tmp_info_state[7] 
                                            + tmp_info_state[1] * tmp_info_state[8]
                                            + tmp_info_state[2] * tmp_info_state[9];
      ukfPrecisionType abs_angle_degree = RadToDeg(acos((ukfPrecisionType) std::abs(scalar_product_m1_m2)));
      
      
      if (rtopSignal >= _rtop_min){
        
        // If the angle between the two fibers is less than 30 and if lambda2(tensor2) < lambda2(tensor1)
        // Then we swap tensor 2 and tensor 1 in the state
        if (abs_angle_degree <= min_angle && tmp_info_state[13] < tmp_info_state[6]) {   
          // Swap     
          SwapState2T(tmp_info_state, p);
        }
        
        // Create the opposite seed
        InverseStateDiffusionPropagator(tmp_info_state, tmp_info_inv_state);
        
        // Update the original directions
        info.start_dir << tmp_info_state[0], tmp_info_state[1], tmp_info_state[2];
        info_inv.start_dir << -tmp_info_state[0], -tmp_info_state[1], -tmp_info_state[2];
        
        // Add the primary seeds to the vector
        info.state = ConvertVector<stdVecState, State>(tmp_info_state);
        info_inv.state = ConvertVector<stdVecState, State>(tmp_info_inv_state);
        seed_infos.push_back(info);
        seed_infos.push_back(info_inv);
        
        // If we have two fibers, we create a second starting seed
        if (abs_angle_degree > min_angle) {
          
          
          // Swap the first and second direction, in the state as well as the covariance
          SwapState2T(tmp_info_state, p);
          
          // Create a crossing seed
          SeedPointInfo info_branching;
          info_branching.point = starting_points[i];
          info_branching.start_dir << tmp_info_state[0], tmp_info_state[1], tmp_info_state[2];
          info_branching.fa = rtop1;
          info_branching.fa2 = rtop2;
          info_branching.trace = rtop;
          info_branching.trace2 = rtopSignal;
          info_branching.covariance.resize(state_dim, state_dim);
          info_branching.covariance = info.covariance;
          
          info_branching.state = ConvertVector<stdVecState, State>(tmp_info_state);       
          seed_infos.push_back(info_branching);
          
          // Create the seed for the opposite direction, keep the other parameters as set for the first direction
          InverseStateDiffusionPropagator(tmp_info_state, tmp_info_inv_state);
          info_branching.state = ConvertVector<stdVecState, State>(tmp_info_inv_state);
          info_branching.start_dir <<  tmp_info_inv_state[0], tmp_info_inv_state[1], tmp_info_inv_state[2];
          seed_infos.push_back(info_branching);          
        }
               
      } // End of if (rtopSignal >= _rtop_min)
      
    } // If not _noddi and not _diffusionPropagator
    else
    {
      tmp_info_state[3] = param[6];     // l1
      tmp_info_state[4] = param[7];     // l2
      tmp_info_inv_state[3] = param[6]; // l1
      tmp_info_inv_state[4] = param[7]; // l2
    }
    // Duplicate/tripple states if we have several tensors.
    if( _num_tensors > 1 && !_noddi && ! _diffusionPropagator)
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
      if( _free_water && ! _diffusionPropagator )
        {
        tmp_info_state.push_back(1);
        tmp_info_inv_state.push_back(1); // add the weight to the state (well was sich rhymt das stiimt)
        }
    }

    if (! _diffusionPropagator) {
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
}

void Tractography::PrintState(State& state)
{
  std::cout<<"State \n";
  std::cout<<"\t m1: "<<state[0]<<" "<<state[1]<<" "<<state[2]<<std::endl;
  std::cout<<"\t l11 .. l14: "<<state[3]<<" "<<state[4]<<" "<<state[5]<<" "<<state[6]<<std::endl;
  std::cout<<"\t m2: "<<state[7]<<" "<<state[8]<<" "<<state[9]<<std::endl;
  std::cout<<"\t l21 .. l24: "<<state[10]<<" "<<state[11]<<" "<<state[12]<<" "<<state[13]<<std::endl;
  std::cout<<"\t w: "<<state[14]<<std::endl;
  std::cout<<" --- "<<std::endl;
}

void Tractography::InverseStateDiffusionPropagator(stdVecState& reference, stdVecState& inverted)
{
  assert(reference.size() == 15);
  assert(inverted.size() == 15);
  for (unsigned int it = 0; it < reference.size(); ++it) {
    if (it <= 2) {
      inverted[it] = -reference[it];
    }
    else {
      inverted[it] = reference[it];
    }
  } 
}

void Tractography::StateToMatrix(State& state, ukfMatrixType& matrix) 
{
  assert(state.size() > 0);
  
  matrix.resize(state.size(), 1);
  
  for (int it = 0; it<state.size(); ++it) {
    matrix(it, 0) = state[it];
  }
  
}

void Tractography::MatrixToState(ukfMatrixType& matrix, State& state) 
{
  assert(matrix.cols() == 1);
  assert(matrix.rows() > 0);
  
  state.resize(matrix.rows());
  
  for (int it = 0; it<matrix.rows(); ++it) {
    state[it] = matrix(it, 0);
  }
  
}

void Tractography::InitLoopUKF(State& state, ukfMatrixType& covariance, ukfVectorType& signal, ukfPrecisionType& dNormMSE)
{

  // INITIALIZATION
  
  // Suppose that the signal is already loaded in memory
  assert(_signal_data);
  
  // For the initialization, we use the model with different parameters, which are hard-coded here.
  // This is the reason why we declare init_model here, and that we do not use _model
  
  ukfPrecisionType l_Qm = 0.005;
  ukfPrecisionType l_Ql = 1500;
  ukfPrecisionType l_Qt = 1500;
  ukfPrecisionType l_Qw = 0.005;
  ukfPrecisionType l_Rs = 0.001;

  
  ukfVectorType weightsOnTensors(2);
    weightsOnTensors[0] = 0.5;
    weightsOnTensors[1] = 0.5;
  
  const bool freeWater = true;  
  ukfPrecisionType D_ISO = 0.003;
  FilterModel* init_model = NULL;
  init_model = new Simple2T_BiExp_FW(l_Qm, l_Ql, l_Qw, l_Rs, weightsOnTensors, freeWater, D_ISO, l_Qt);
  
  if (init_model == NULL) {
    std::cout<<"NEW operator failed to instentiate a Simple2T_BiExp_FW model during initialization. Exiting the program"<<std::endl;
    exit(1);
  }
  
  init_model->set_signal_data(_signal_data);
  init_model->set_signal_dim(_signal_data->GetSignalDimension() * 2);
  
  UnscentedKalmanFilter* ukf = NULL;
  ukf = new UnscentedKalmanFilter(init_model);
  
  if (init_model == NULL) {
    std::cout<<"NEW operator failed to instentiate an UnscentedKalmanFilter during initialization. Exiting the program"<<std::endl;
    exit(1);
  }
  

  ukfPrecisionType dNormMSE_min;
  ukfPrecisionType er;
  ukfPrecisionType er_min;
  
  
  State state_new( init_model->state_dim() );
  State state_min( init_model->state_dim() );

  ukfMatrixType state_matrix(init_model->state_dim(), 1);
  
  ukfMatrixType covariance_new( init_model->state_dim(), init_model->state_dim() );
  covariance_new.setConstant(ukfZero);  
  
  ukfMatrixType covariance_min( init_model->state_dim(), init_model->state_dim() );
  
  state_min = state;
  covariance_min = covariance;
  

  // END OF INITIALIZATION
  
  er = 0.1; 
  state_min = state;
  covariance_min = covariance;
  dNormMSE_min = dNormMSE;
  
  int max_iter = 100;
  //int max_iter = 3;
  int jj  = 1;
  er_min = er;
  
  while (er > 0.003) {

    ukf->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE); 
    
    // 1. F function
    StateToMatrix(state_new, state_matrix);
    init_model->F(state_matrix);
    
    // 2. Compute error
    er = dNormMSE;
    
    // 3. Update state
    MatrixToState(state_matrix, state_new);    
    state = state_new;
    covariance = covariance_new;
    
    // 4. Look at the error
    if (er < er_min) {
      state_min = state;
      covariance_min = covariance;
      er_min = er;
      dNormMSE_min = dNormMSE;
    }
    
    if (jj > max_iter) {
      // std::cout<<"Maximum number of iteration reached"<<std::endl;
      break;
    }
       
    ++jj;
  }
  
  state = state_min;
  covariance = covariance_min;
  dNormMSE = dNormMSE_min;
  

  delete init_model;
  delete ukf;

}

void Tractography::computeRTOPfromSignal(ukfPrecisionType& rtopSignal, ukfVectorType& signal)
{
  assert(signal.size() > 0);
  
  rtopSignal = 0.0;
  
  // The RTOP is the sum of the signal
  // We use signal.size()/2 because the first half of the signal is identical
  // to the second half.
  
  for (int i=0; i < signal.size()/2; ++i) {
    rtopSignal += signal[i];
    
    if (signal[i] < 0) {
      std::cout<<"Negative signal found when computing the RTOP from the signal, value : "<<signal[i]<<std::endl;
    }
    
  }  
}


void Tractography::computeRTOPfromState(stdVecState& state, ukfPrecisionType& rtop, ukfPrecisionType& rtop1, ukfPrecisionType& rtop2)
{
  // Control input: state should have 15 rows
  assert(state.size() == 15);
  
  ukfPrecisionType l11 = state[3]*1e-6;
  ukfPrecisionType l12 = state[4]*1e-6;
  ukfPrecisionType l13 = state[5]*1e-6;
  ukfPrecisionType l14 = state[6]*1e-6;
  
  ukfPrecisionType l21 = state[10]*1e-6;
  ukfPrecisionType l22 = state[11]*1e-6;
  ukfPrecisionType l23 = state[12]*1e-6;
  ukfPrecisionType l24 = state[13]*1e-6;
  
  ukfPrecisionType w = state[14];
  
  ukfPrecisionType det_l1 = l11 * l12;
  ukfPrecisionType det_t1 = l13 * l14;
  ukfPrecisionType det_l2 = l21 * l22;
  ukfPrecisionType det_t2 = l23 * l24;
  
  // !!! D_ISO value hardcoded...
  ukfPrecisionType det_fw = 0.003*0.003*0.003;
  
  // !!! 0.7 and 0.3 tensor weights are hardcoded...
  rtop1 = std::pow(UKF_PI, 3/2) * w * (0.7/std::sqrt(det_l1) + 0.3/std::sqrt(det_t1));
  rtop2 = std::pow(UKF_PI, 3/2) * w * (0.7/std::sqrt(det_l2) + 0.3/std::sqrt(det_t2));
  rtop  = rtop1 + rtop2 + 2*std::pow(UKF_PI, 3/2) * (1-w)/std::sqrt(det_fw);
  
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
  if (static_cast<int>(primary_seed_infos.size()) < 1) {
    std::cout<<"The number of seeds after the initialization step is less than one (1). Cannot perform tractography. Exiting"<<std::endl;
    exit(1);
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
    std::cout << "Tracing " << primary_seed_infos.size() << " primary fibers:" << std::endl;
    raw_primary.resize(primary_seed_infos.size() );
    // Timer timer ;
    WorkDistribution work_distribution = GenerateWorkDistribution(num_of_threads,
                                                                  static_cast<int>(primary_seed_infos.size() ) );

    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( num_of_threads );
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
      threader->SetMultipleMethod(i, ThreadCallback, &str);
      }
    threader->SetGlobalDefaultNumberOfThreads(num_of_threads);
    threader->MultipleMethodExecute();

    // std::cout << "Time cost: " << timer.elapsed() << std::endl << std::endl ;

    // Unpack the branch seeds and their affiliation
    int num_branch_seeds = 0;
    for( int i = 0; i < num_of_threads; i++ )
      {
      num_branch_seeds += static_cast<int>(str.branching_seed_info_vec->at(i).size() );
      }

    std::cout << "branch_seeds size: " << num_branch_seeds << std::endl;
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
    std::cout << "Tracing " << branch_seed_infos.size() << " branches:" << std::endl;

    raw_branch.resize(branch_seed_infos.size() );
    // Timer timer ;
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

    itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( num_of_threads );
    for( int i = 0; i < num_of_threads; i++ )
      {
      threader->SetMultipleMethod(i, ThreadCallback, &str);
      }
    threader->MultipleMethodExecute();

    // std::cout << "Time cost: " << timer.elapsed() << std::endl << std::endl ;
    }

  std::vector<UKFFiber> fibers;
  PostProcessFibers(raw_primary, raw_branch, branch_seed_affiliation, _branches_only, fibers);
  std::cout << "fiber size after postprocessibers: " << fibers.size() << std::endl;
  if ( fibers.size()  == 0 )
  {
    return EXIT_FAILURE;
  }

  // Write the fiber data to the output vtk file.
  VtkWriter writer(_signal_data, _filter_model_type, _record_tensors);
  // possibly write binary VTK file.
  writer.SetWriteBinary(this->_writeBinary);
  writer.SetWriteCompressed(this->_writeCompressed);

  writer.set_transform_position(_transform_position);
  const int writeStatus = writer.Write(_output_file, _output_file_with_second_tensor, fibers, _record_state, _store_glyphs, _noddi, _diffusionPropagator);
  // Clear up the kalman filters
  for( size_t i = 0; i < _ukf.size(); i++ )
    {
    delete _ukf[i];
    }
  _ukf.clear();
  return writeStatus;
}

void Tractography::createProtocol(const ukfVectorType& _b_values,
                    ukfVectorType& _gradientStrength, ukfVectorType& _pulseSeparation)
{
  std::vector<double> Bunique, tmpG;
  ukfPrecisionType Bmax = 0;
  ukfPrecisionType tmp, Gmax, GAMMA;

  _gradientStrength.resize(_b_values.size());
  _pulseSeparation.resize(_b_values.size());

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
    _pulseSeparation[i] = tmp;
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
        _gradientStrength[j] = tmpG[i];
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

  std::cout << std::endl;

  assert(ret.size() == s.size() );

  // Build B matrix.
  const int signal_dim = _signal_data->GetSignalDimension();
 int B_size = 0;
  
  // We define the limits for the diffusionPropagator
  const int HIGH_BVAL = 1500;
  const int LOW_BVAL = 800;
  
  // When we use the diffusion propagator, we would like to take into account only the gradients with a BValue between 800 and 1500 for the initialisation
  if (_diffusionPropagator) {
    
    // We count the number of gradient directions with LOW_BVAL < bval <= HIGH_BVAL
    for (int i = 0; i < signal_dim * 2; ++i) {
      if (b[i] > LOW_BVAL && b[i] <= HIGH_BVAL) {
        B_size++;
      }
    }
    
  } 
  else {
    B_size = signal_dim * 2;
  }
  
  std::cout<<"Using "<<B_size/2<<" out of "<<signal_dim<<" b-values during initialization, when estimating a single tensor"<<std::endl;
  
/**
  * A special type for holding 6 elements of tensor for each signal
  *  Only used in tractography.cc
  */
  typedef Eigen::Matrix<ukfPrecisionType, Eigen::Dynamic, 6>  BMatrixType;
  BMatrixType B(B_size, 6); //HACK: Eigen::Matrix<ukfPrecisionType, DYNAMIC, 6> ??

  int local_counter = 0;
  for( int i = 0; i < signal_dim * 2; ++i )
    {
      if ((b[i] <= LOW_BVAL || b[i] > HIGH_BVAL) && (_diffusionPropagator)) {
        continue;
      }
    const vec3_t & g = u[i];
    B(local_counter, 0) = (-b[i]) * (g[0] * g[0]);
    B(local_counter, 1) = (-b[i]) * (2.0 * g[0] * g[1]);
    B(local_counter, 2) = (-b[i]) * (2.0 * g[0] * g[2]);
    B(local_counter, 3) = (-b[i]) * (g[1] * g[1]);
    B(local_counter, 4) = (-b[i]) * (2.0 * g[1] * g[2]);
    B(local_counter, 5) = (-b[i]) * (g[2] * g[2]);
    
    local_counter++;
    }
  // Use QR decomposition to find the matrix representation of the tensor at the
  // seed point of the fiber. Raplacement of the gmm function gmm::least_squares_cg(..)
  Eigen::HouseholderQR<BMatrixType> qr(B);

  // Temporary variables.
  mat33_t D;

  std::cout << "Estimating seed tensors:" << std::endl;
  // Timer timer ;
  // boost::progress_display disp(static_cast<unsigned long>(ret.size())) ;
  // Unpack data
  for( stdEigVec_t::size_type i = 0; i < s.size(); ++i )
    {
    // We create a temporary vector to store the signal when the log function is applied
    ukfVectorType log_s;
    log_s.resize(s[i].size());
    
    for( unsigned int j = 0; j < s[i].size(); ++j )
      {
      if( s[i][j] <= 0 )
        {
        s[i][j] = 10e-8;
        }
      //s[i][j] = log(s[i][j]);
      log_s[j] = log(s[i][j]);
      }

    // The six tensor components.
    //TODO: this could be fixed size
    // ukfVectorType d = qr.solve(s[i]);
    ukfVectorType d = qr.solve(log_s);

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

    // ++disp ; for boost progress bar
    }

  // std::cout << "Time cost: " << timer.elapsed() << std::endl << std::endl ;
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
    const bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1;

    state_tmp.col(0) = state;
    _model->H(state_tmp, signal_tmp);

    const ukfPrecisionType ga = s2ga(signal_tmp);
    const bool   in_csf = ga < _ga_min || fa < _fa_min;
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
  ukfMatrixType signal_tmp(_model->signal_dim(), 1); // estimated signal from the state
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
    const bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1; // is this 0.1 correct? yes

    // after here state does not change until next step.

    state_tmp.col(0) = state;

    _model->H(state_tmp, signal_tmp); // signal_tmp is written, but only used to calculate ga

    const ukfPrecisionType ga = s2ga(signal_tmp);
    bool in_csf = (_noddi) ? ( ga < _ga_min ) : (ga < _ga_min || fa < _fa_min);

    bool dNormMSE_too_high = false;
    bool negative_free_water = false;
    
    if (_diffusionPropagator) {
      ukfPrecisionType rtopSignal = trace2; // rtopSignal is stored in trace2
      in_csf = rtopSignal <  _rtop_min;
      dNormMSE_too_high = dNormMSE > _maxNMSE;
      negative_free_water = state[14] < 0.0; 
    }
    const bool is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain
        || in_csf
        || stepnr > _max_length  // Stop when the fiber is too long
        || is_curving
        || dNormMSE_too_high
        || negative_free_water)
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
    if( is_branching && ! _diffusionPropagator )
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
      
      else if (is_branching && _diffusionPropagator) {
        ukfPrecisionType scalar_product_m1_m2 = m1.dot(m2); // m1 and m2 are unit vectors
        ukfPrecisionType angle = RadToDeg(acos(scalar_product_m1_m2));
        
        ukfPrecisionType rtop_signal = trace2;
        ukfPrecisionType rtop2 = fa2;

        // If we have a crossing fiber (angle greater than 25 degree)
        if (rtop_signal >= _rtop_min && rtop2 >= _rtop_min && angle >= 25) {
          
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
          local_seed.fa = rtop2;

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
    const bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1; // x is the seed point
    state_tmp.col(0) = state;

    _model->H(state_tmp, signal_tmp);

    const ukfPrecisionType ga = s2ga(signal_tmp);
    bool in_csf;
    if(_noddi)
      in_csf = ga < _ga_min;
    else
      in_csf = ga < _ga_min || fa < _fa_min;
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

void Tractography::LoopUKF(const int thread_id,
                           State& state, 
                           ukfMatrixType& covariance, 
                           ukfVectorType& signal, 
                           State& state_new, 
                           ukfMatrixType& covariance_new, 
                           ukfPrecisionType& dNormMSE
                           )
{  

  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);     
  state = state_new;
  covariance = covariance_new;
  
  ukfPrecisionType er_org = 0.0;
  ukfPrecisionType er = 0.0;
  
  er_org = dNormMSE;
  er = er_org;
  
  State state_prev = state;
  
  int max_iter = _maxUKFIterations;
  
  for (int jj = 0; jj < max_iter; ++jj) {
    _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);     
    state = state_new;
    //covariance = covariance_new;
    
    er_org = er;
    er = dNormMSE;
      
    if (er_org - er < 0.001) {
      break;
    }
    
    state_prev = state;
  }
  
  state = state_prev;
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

  if (_diffusionPropagator) {
    LoopUKF(thread_id, state, covariance, signal, state_new, covariance_new, dNormMSE);
  }
  else {
    _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);
    state = state_new;
    covariance = covariance_new;
  }

  const vec3_t old_dir = m1;   // Direction in last step
  ukfPrecisionType fa_tensor_1;
  ukfPrecisionType fa_tensor_2;
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
  else if (_diffusionPropagator) {
    // lxx are not used
    vec3_t l11, l12, l21, l22;
    // we do not use l1 and l2 in the case of the diffusionPropagator model, so we set the values to 0.
    l1 << 0.0, 0.0, 0.0;
    l2 << 0.0, 0.0, 0.0;
    
    // DEBUG
    // Apply F function to normalize directions and bound parameters
    /*ukfMatrixType state_matrix; StateToMatrix(state, state_matrix);
    _model->F(state_matrix);
    MatrixToState(state_matrix, state);*/
    
    _model->State2Tensor2T(state, old_dir, m1, l11, m2, l21);
    
    /*vec3_t dir; 
    initNormalized(dir, state[0], state[1], state[2]);
    state[0] = dir[0]; state[1] = dir[1]; state[2] = dir[2];
    initNormalized(dir, state[7], state[8], state[9]);
    state[7] = dir[0]; state[8] = dir[1]; state[9] = dir[2];*/
    
    ukfPrecisionType rtop1, rtop2, rtopModel, rtopSignal;
    stdVecState local_state = ConvertVector<State, stdVecState>(state);
    
    computeRTOPfromState(local_state, rtopModel, rtop1, rtop2);
    computeRTOPfromSignal(rtopSignal, signal);
    
    fa = rtop1;
    fa2 = rtop2;
    trace = rtopModel;
    trace2 = rtopSignal;
    
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
    else if (_diffusionPropagator) 
      {
        /*vec3_t avg_m;
        initNormalized(avg_m, m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2]);
        state[0] = state[7] = avg_m[0];
        state[1] = state[8] = avg_m[1];
        state[2] = state[9] = avg_m[2];
        m1 = avg_m;
        m2 = avg_m;*/
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
  // exit(1);
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

void Tractography::SwapState2T(stdVecState& state, ukfMatrixType& covariance)
{
  State tmp_state = ConvertVector<stdVecState, State>(state);
  SwapState2T(tmp_state, covariance);
  state = ConvertVector<State, stdVecState>(tmp_state);
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

    
  if( (_record_trace || _record_kappa) && !_recordRTOP)
    {
    fiber.trace.push_back(2*(atan(1/trace)/3.14));
    if( _num_tensors >= 2 )
      {
      fiber.trace2.push_back(2*(atan(1/trace2)/3.14));
      }
    }
    
  if (_recordRTOP){
    fiber.trace.push_back(trace);
    fiber.trace2.push_back(trace2);
  }

  if( _record_fa || _record_Vic || _recordRTOP)
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
        exit(1);
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
        exit(1);
        }
      }
    fiber.free_water.push_back(fw);
    }

  // Record the state
  if( (state.size() == 5 || state.size() == 10 || state.size() == 15) && ! _diffusionPropagator)   // i.e. simple model
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
  else if (_diffusionPropagator) {
    State store_state(state);
    vec3_t dir;
    
    // normalize m1
    initNormalized(dir, store_state[0], store_state[1], store_state[2]);
    store_state[0] = dir[0];
    store_state[1] = dir[1];
    store_state[2] = dir[2];    
    
    // normalize m2
    initNormalized(dir, store_state[7], store_state[8], store_state[9]);
    store_state[7] = dir[0];
    store_state[8] = dir[1];
    store_state[9] = dir[2];   
    
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

  if( _record_trace || _record_kappa || _recordRTOP)
    {
    fiber.trace.reserve(fiber_size);
    if( _num_tensors >= 2 )
      {
      fiber.trace2.reserve(fiber_size);
      }
    }

  if( _record_fa || _record_Vic || _recordRTOP)
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

  // Calculate the Normalized Mean Square Error between
  // the signal and the estimated parameters
itk::SingleValuedCostFunction::MeasureType itk::DiffusionPropagatorCostFunction::GetValue(const ParametersType &parameters) const
{
  MeasureType residual = 0.0;
    
  // TODO: check parameters size
  
  // Convert the parameter to the ukfMtarixType
  ukfMatrixType localState(this->GetNumberOfParameters(), 1);
  for (unsigned int it=0; it<this->GetNumberOfParameters(); ++it) {
    localState(it, 0) = parameters[it];
  }
  // w for freeWater
  if (localState(14, 0) < 0.0) {
    localState(14, 0) = 0.0;
  }    
  // Estimate the signal
  ukfMatrixType estimatedSignal(this->GetNumberOfValues(), 1);

  //_model->F(localState);
  _model->H(localState, estimatedSignal);
  
  // Compute the error between the estimated signal and 
  // the acquired one
  ukfPrecisionType err = 0.0;
  this->computeError(estimatedSignal, _signal, err);
    
  // Return the result    
  residual = err;    
  return residual;
}
  
void itk::DiffusionPropagatorCostFunction::GetDerivative(const ParametersType &parameters, DerivativeType & derivative ) const 
{

  // We use numerical derivative
  // slope = [f(x+h) - f(x-h)] / (2h)
    
  ParametersType p_h(this->GetNumberOfParameters());    // for f(x+h)
  ParametersType p_hh(this->GetNumberOfParameters());   // for f(x-h)
  
  // The size of the derivative is not set by default,
  // so we have to do it manually
  derivative.SetSize(this->GetNumberOfParameters());

  // Set parameters
  for (unsigned int it=0; it<this->GetNumberOfParameters(); ++it) {
    p_h[it] = parameters[it];
    p_hh[it] = parameters[it];  
  }
  
  // Calculate derivative for each parameter
  for (unsigned int it=0; it<this->GetNumberOfParameters(); ++it) {
  
    // Optimal h is sqrt(epsilon machine) * x
    double h = std::sqrt(2.3e-16)*std::max(std::abs(parameters[it]), 1.0);
    
    // Volatile, otherwise compiler will optimize the value for dx
    volatile double xph = parameters[it] + h;
    
    // For taking into account the rounding error
    double dx = xph - parameters[it];
    
    // Compute the slope
    p_h[it] = xph;
    //p_hh[it] = parameters[it] - h;  
    derivative[it] = (this->GetValue(p_h) - this->GetValue(p_hh))/(dx);
    
    // Set parameters back for next iteration
    p_h[it] = parameters[it];
    p_hh[it] = parameters[it];
  }
    
}


#include "itkLBFGSBOptimizer.h"

void Tractography::NonLinearLeastSquareOptimization(State& state, ukfVectorType& signal, FilterModel *model) 
{
  typedef itk::LBFGSBOptimizer                  OptimizerType;
  typedef itk::DiffusionPropagatorCostFunction  CostType;

  CostType::Pointer cost = CostType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();

  cost->SetNumberOfParameters(state.size());
  cost->SetNumberOfValues(signal.size());
  cost->SetSignalValues(signal);
  cost->SetModel(model);

  optimizer->SetCostFunction(cost);
  
  CostType::ParametersType p(cost->GetNumberOfParameters());
  // Fill p
  for (int it=0; it<state.size(); ++it) {
     p[it] = state[it];
  }  
  optimizer->SetInitialPosition(p);
  optimizer->SetProjectedGradientTolerance(1e-10);
  optimizer->SetMaximumNumberOfIterations(500);
  optimizer->SetMaximumNumberOfEvaluations(500);
  optimizer->SetMaximumNumberOfCorrections(5); // ? Meaning ?
  optimizer->SetCostFunctionConvergenceFactor(1e2); // Precision of the solution 1e1 = very precise, 1e12 = coarse estimate
  optimizer->SetTrace(false); // Print debug info
  
  // Set bounds
  
  OptimizerType::BoundSelectionType boundSelect(cost->GetNumberOfParameters());
  OptimizerType::BoundValueType upperBound(cost->GetNumberOfParameters());
  OptimizerType::BoundValueType lowerBound(cost->GetNumberOfParameters());  
  
  boundSelect.Fill(2); // BOTHBOUNDED = 2
  lowerBound.Fill(0.0);
  upperBound.Fill(3000.0);
  
  // Lower bound
  // First bi-exponential parameters
  lowerBound[0] = lowerBound[1] = lowerBound[2] = -1.0;
  lowerBound[3] = lowerBound[4] = 1.0;
  lowerBound[5] = lowerBound[6] = 0.1;

  // Second bi-exponential
  lowerBound[7] = lowerBound[8] = lowerBound[9] = -1.0;
  lowerBound[10] = lowerBound[11] = 1.0;  
  lowerBound[12] = lowerBound[13] = 0.1;   
  lowerBound[14] = 0.0; // free water between 0 and 1
   
  // Upper bound
  // First bi-exponential
  upperBound[0] = upperBound[1] = upperBound[2] = 1.0;
  upperBound[3] = upperBound[4] = upperBound[5] = upperBound[6] = 3000.0;
  
  // Second bi-exponential
  upperBound[7] = upperBound[8] = upperBound[9] = 1.0;
  upperBound[10] = upperBound[11] = upperBound[12] = upperBound[13] = 3000.0;
  upperBound[14] = 1.0;
  
  optimizer->SetBoundSelection(boundSelect);
  optimizer->SetUpperBound(upperBound);
  optimizer->SetLowerBound(lowerBound);  
  
  
  optimizer->StartOptimization();

  // std::cout<<"Current Iteration: "<<optimizer->GetCurrentIteration()<<std::endl; 
  // std::cout<<"Value: "<<optimizer->GetValue()<<std::endl;
  p = optimizer->GetCurrentPosition();
  // Write back the state
  for (int it=0; it<state.size(); ++it) {
    state[it] = p[it];
  }

}

