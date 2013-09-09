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

#include <vnl/algo/vnl_determinant.h>
#include <fstream>
#include <iostream>

// TODO implement this switch
#include "config.h"

Tractography::Tractography(FilterModel *model, model_type filter_model_type,

                           const std::string& output_file, const std::string & output_file_with_second_tensor,
                           const bool record_fa, const bool record_nmse, const bool record_trace,
                           const bool record_state,
                           const bool record_cov, const bool record_free_water,  const bool record_tensors,
                           const bool transform_position, const bool store_glyphs, const bool branchesOnly,

                           const double fa_min, const double ga_min, const double seedFALimit,
                           const int num_tensors, const int seeds_per_voxel,
                           const double minBranchingAngle, const double maxBranchingAngle,
                           const bool is_full_model, const bool free_water,
                           const double stepLength, const double maxHalfFiberLength,
                           const std::vector<int>& labels,

                           double p0, double sigma_signal, double sigma_mask,
                           double min_radius, double full_brain_ga_min,

                           const int num_threads
                           ) :
  _ukf(0, NULL), _model(model), _filter_model_type(filter_model_type),

  _output_file(output_file), _output_file_with_second_tensor(output_file_with_second_tensor),
  _record_fa(record_fa), _record_nmse(record_nmse), _record_trace(record_trace), _record_state(record_state),
  _record_cov(record_cov), _record_free_water(record_free_water), _record_tensors(record_tensors),
  _transform_position(transform_position), _store_glyphs(store_glyphs), _branches_only(branchesOnly),

  _p0(p0), _sigma_signal(sigma_signal), _sigma_mask(sigma_mask), _min_radius(min_radius), _full_brain_ga_min(
    full_brain_ga_min),
  _max_length(static_cast<int>(std::ceil(maxHalfFiberLength / stepLength) ) ), _full_brain(false),

  _fa_min(fa_min), _ga_min(ga_min), _seedFALimit(seedFALimit),
  _num_tensors(num_tensors), _seeds_per_voxel(seeds_per_voxel),
  _cos_theta_min(minBranchingAngle), _cos_theta_max(maxBranchingAngle),
  _is_full_model(is_full_model), _free_water(free_water),
  _stepLength(stepLength),
  _labels(labels),
  _writeBinary(false),
  _num_threads(num_threads)
{
  if( _cos_theta_max != 0.0 && _cos_theta_max <= _cos_theta_min )
    {
    std::cout << "Maximum branching angle must be greater than " << minBranchingAngle << " degrees." << std::endl;
    exit(1);
    }

  if( _num_tensors < 1 || _num_tensors > 3 )
    {
    std::cout << "Only one, two or three tensors are supported." << std::endl;
    exit(1);
    }

  _cos_theta_max = cos(_cos_theta_max * M_PI / 180.0);
  _cos_theta_min = cos(_cos_theta_min * M_PI / 180.0);

  // Double check branching.
  _is_branching = _num_tensors > 1 && _cos_theta_max < 1.0;  // The branching is enabled when the maximum branching
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
    }
  else if( _num_tensors == 2 )    // 2 TENSOR CASE ////////////////////////////////
    {
    if( !is_full_model && free_water ) // simple model with free water
      {
      _nPosFreeWater = 10;
      }
    else if( is_full_model && free_water ) // full model with free water
      {
      _nPosFreeWater = 12;
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

  std::vector<vec_t> seeds;
  assert(_labels.size() > 0);
  if( !_full_brain )
    {
    _signal_data->GetSeeds(_labels, seeds);
    }
  else
    {
    // Iterate through all brain voxels and take those as seeds voxels.
    const vec_t dim = _signal_data->dim();
    for( int x = 0; x < dim._[0]; ++x )
      {
      for( int y = 0; y < dim._[1]; ++y )
        {
        for( int z = 0; z < dim._[2]; ++z )
          {
          vec_t pos = make_vec(x, y, z);
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
  std::vector<vec_t> rand_dirs;

  if( seeds.size() == 1 && _seeds_per_voxel == 1 )   // if there is only one seed don't use offset so fibers can be
                                                     // compared
    {
    rand_dirs.push_back(make_vec(0, 0, 0) );   // in the test cases.
    }
  else
    {
    for( int i = 0; i < _seeds_per_voxel; ++i )
      {
      vec_t dir = make_vec(static_cast<double>( (rand() % 10001) - 5000),
                           static_cast<double>( (rand() % 10001) - 5000),
                           static_cast<double>( (rand() % 10001) - 5000) );

      // CB: those directions are to compare against the matlab output
      // dir._[2] = 0.439598093988175;
      // dir._[1] = 0.236539281163321;
      // dir._[0] = 0.028331682419209;

      dir /= norm(dir);
      dir *= 0.5;

      rand_dirs.push_back(dir);
      }
    }
  // Calculate all starting points.
  std::vector<vec_t>                starting_points;
  std::vector<std::vector<double> > signal_values;
  std::vector<double>               signal;
  signal.resize(signal_dim * 2);

  std::vector<vec_t>::const_iterator cit;
  std::vector<vec_t>::iterator       jt;
  int                                num_less_than_zero = 0;
  int                                num_invalid = 0;
  int                                num_ga_too_low = 0;

  int tmp_counter = 1;
  for( cit = seeds.begin(); cit != seeds.end(); ++cit )
    {
    for( jt = rand_dirs.begin(); jt != rand_dirs.end(); ++jt )
      {
      vec_t point = *cit + *jt;

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

        if( isnan(signal[k]) || isinf(signal[k]) )
          {
          keep = false;
          ++num_invalid;
          break;
          }

        vnl_vector_ref<double> signal_vnl(signal.size(), &signal.front() );

        // If we do full brain tractography we only want seed voxels where the
        // GA is bigger than 0.18.
        vnl_matrix<double> signal_tmp(signal_dim * 2, 1);
        signal_tmp.set_column(0, signal_vnl);
        if( _full_brain && s2ga(signal_tmp) < _full_brain_ga_min )
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
  std::vector<std::vector<double> > starting_params(starting_points.size() );

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
    const std::vector<double>& param = starting_params[i];

    assert(param.size() == 9);

    // Filter out seeds whose FA is too low.
    double fa = l2fa(param[6], param[7], param[8]);
    double trace = param[6] + param[7] + param[8];
    double fa2 = -1;
    double trace2 = -1;

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

    // Create seed info for both directions;
    SeedPointInfo info;
    SeedPointInfo info_inv;

    info.point = starting_points[i];
    info.start_dir = make_vec(param[0], param[1], param[2]);
    info.fa = fa;
    info.fa2 = fa2;
    info.trace = trace;
    info.trace2 = trace2;
    info_inv.point = starting_points[i];
    info_inv.start_dir = make_vec(-param[0], -param[1], -param[2]);
    info_inv.fa = fa;
    info_inv.fa2 = fa2;
    info_inv.trace = trace;
    info_inv.trace2 = trace2;

    if( _is_full_model )
      {
      info.state.resize(6);
      info_inv.state.resize(6);
      info.state[0] = param[3];     // Theta
      info.state[1] = param[4];     // Phi
      info.state[2] = param[5];     // Psi
      info.state[5] = param[8];     // l3
      info_inv.state[0] = param[3]; // Theta
      info_inv.state[1] = param[4]; // Phi
      // Switch psi angle.
      // Careful here since M_PI is not standard c++.
      info_inv.state[2] = (param[5] < 0.0 ? param[5] + M_PI : param[5] - M_PI);
      info_inv.state[5] = param[8]; // l3

      }
    else     // i.e. simple model
      { // Starting direction.
      info.state.resize(5);
      info_inv.state.resize(5);
      info.state[0] = info.start_dir._[0];
      info.state[1] = info.start_dir._[1];
      info.state[2] = info.start_dir._[2];
      info_inv.state[0] = info_inv.start_dir._[0];
      info_inv.state[1] = info_inv.start_dir._[1];
      info_inv.state[2] = info_inv.start_dir._[2];
      }

    info.state[3] = param[6];     // l1
    info.state[4] = param[7];     // l2
    info_inv.state[3] = param[6]; // l1
    info_inv.state[4] = param[7]; // l2

    // Duplicate/tripple states if we have several tensors.
    if( _num_tensors > 1 )
      {
      size_t size = info.state.size();
      for( size_t j = 0; j < size; ++j )
        {
        info.state.push_back(info.state[j]);
        info_inv.state.push_back(info.state[j]);
        }
      if( _num_tensors > 2 )
        {
        for( size_t j = 0; j < size; ++j )
          {
          info.state.push_back(info.state[j]);
          info_inv.state.push_back(info.state[j]);
          }
        }
      }

    if( _free_water )
      {
      info.state.push_back(1);
      info_inv.state.push_back(1); // add the weight to the state (well was sich rhymt das stiimt)
      }

    int state_dim = info.state.size();

    info.covariance.set_size(state_dim, state_dim);
    info_inv.covariance.set_size(state_dim, state_dim);

    // make sure covariances are really empty
    info.covariance.fill(0);
    info_inv.covariance.fill(0);
    for( int local_i = 0; local_i < state_dim; ++local_i )
      {
      info.covariance(local_i, local_i) = _p0;
      info_inv.covariance(local_i, local_i) = _p0;
      }

    seed_infos.push_back(info);
    seed_infos.push_back(info_inv);   // NOTE that the seed in reverse direction is put directly after the seed in
                                      // original direction
    }
}

void Tractography::Run()
{
  assert(_signal_data);   // The _signal_data is initialized in Tractography::LoadFiles(),
  // Thus Run() must be invoked after LoadFiles()
  // Initialize and prepare seeds.

  std::vector<SeedPointInfo>            primary_seed_infos;
  std::vector<SeedPointInfo>            branch_seed_infos;       // The info of branching seeds
  std::vector<BranchingSeedAffiliation> branch_seed_affiliation; // Which fiber originated from the main seeds is this
                                                                 // branch attached

  Init(primary_seed_infos);

  const int num_of_threads = std::min(_num_threads, static_cast<int>(primary_seed_infos.size() ) );
  // const int num_of_threads = 8;

  assert(num_of_threads > 0);
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

  // Write the fiber data to the output vtk file.
  VtkWriter writer(_signal_data, _filter_model_type, _record_tensors);
  // possibly write binary VTK file.
  writer.SetWriteBinary(this->_writeBinary);

  writer.set_transform_position(_transform_position);
  writer.Write(_output_file, _output_file_with_second_tensor, fibers, _record_state, _store_glyphs);
  // Clear up the kalman filters
  for( size_t i = 0; i < _ukf.size(); i++ )
    {
    delete _ukf[i];
    }
  _ukf.clear();
}

void Tractography::UnpackTensor(const std::vector<double>& b,
                                const std::vector<vec_t>& u,            // u - directions
                                std::vector<std::vector<double> >& s,   // s = signal values
                                std::vector<std::vector<double> >& ret) // starting params [i][j] : i - signal number; j
                                                                        // - param
{

  // DEBUGGING
  // std::cout << "b's: ";
  // for (int i=0; i<b.size();++i) {
  //   std::cout << b[i] << ", ";
  // }

  std::cout << std::endl;

  int signal_dim = _signal_data->GetSignalDimension();

  assert(ret.size() == s.size() );

  // Build B matrix.
  vnl_matrix<double> B(signal_dim * 2, 6);
  for( int i = 0; i < signal_dim * 2; ++i )
    {
    const vec_t& g = u[i];
    B(i, 0) = (-b[i]) * (g._[0] * g._[0]);
    B(i, 1) = (-b[i]) * (2.0 * g._[0] * g._[1]);
    B(i, 2) = (-b[i]) * (2.0 * g._[0] * g._[2]);
    B(i, 3) = (-b[i]) * (g._[1] * g._[1]);
    B(i, 4) = (-b[i]) * (2.0 * g._[1] * g._[2]);
    B(i, 5) = (-b[i]) * (g._[2] * g._[2]);

    }

  // The six tensor components.
  vnl_vector<double> d(6);

  // Temporary variables.
  vnl_matrix<double> D(3, 3);
  vnl_matrix<double> Q(3, 3);
  vnl_matrix<double> QT(3, 3);
  vnl_vector<double> sigma(3);
  double             theta, phi, psi;

  std::cout << "Estimating seed tensors:" << std::endl;
  // Timer timer ;
  // boost::progress_display disp(static_cast<unsigned long>(ret.size())) ;
  // Unpack data
  for( size_t i = 0; i < s.size(); ++i )
    {
    for( size_t j = 0; j < s[i].size(); ++j )
      {
      if( s[i][j] <= 0 )
        {
        s[i][j] = 10e-8;

        }

      s[i][j] = log(s[i][j]);
      }

    // Use QR decomposition to find the matrix representation of the tensor at the
    // seed point of the fiber. Raplacement of the gmm function gmm::least_squares_cg(..)
    vnl_vector_ref<double> s_i_vnl(s[i].size(), &s[i].front() );

    vnl_qr<double> QR(B);

    d = QR.solve(s_i_vnl);

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

    vnl_svd<double> svd_decomp(D);
    Q = svd_decomp.U();
    // QT = (svd_decomp.V()).transpose(); // NOTE: QT is never used

    sigma = svd_decomp.W().diagonal(); // diagonal() returns elements of a diag matrix as a vector.

    assert(sigma[0] >= sigma[1] && sigma[1] >= sigma[2]);
    if( vnl_determinant(Q) < 0 )
      {
      Q = Q * (-1.0);
      }
    assert(vnl_determinant(Q) > 0);

    // Extract the three Euler Angles from the rotation matrix.
    theta = acos(Q(2, 2) );
    double epsilon = 1.0e-10;
    if( fabs(theta) > epsilon )
      {
      phi = atan2(Q(1, 2), Q(0, 2) );
      psi = atan2(Q(2, 1), -Q(2, 0) );
      }
    else
      {
      phi = atan2(-Q(0, 1), Q(1, 1) );
      psi = 0.0;
      }

    ret[i].resize(9);
    ret[i][0] = Q(0, 0);
    ret[i][1] = Q(1, 0);
    ret[i][2] = Q(2, 0);
    ret[i][3] = theta;
    ret[i][4] = phi;
    ret[i][5] = psi;
    sigma = sigma * 1.0e6; // NOTICE this scaling of eigenvalues. The values are scaled back in diffusion_euler()
    ret[i][6] = sigma[0];
    ret[i][7] = sigma[1];
    ret[i][8] = sigma[2];

    // ++disp ; for boost progress bar
    }

  // std::cout << "Time cost: " << timer.elapsed() << std::endl << std::endl ;
}

void Tractography::Follow3T(const int thread_id,
                            const size_t seed_index,
                            const SeedPointInfo& seed,
                            UKFFiber& fiber,
                            bool is_branching,
                            std::vector<SeedPointInfo>& branching_seeds,
                            std::vector<BranchingSeedAffiliation>& branching_seed_affiliation)
{

  assert(_model->signal_dim() == _signal_data->GetSignalDimension() * 2);

  // Unpack the seed information.
  vec_t              x = seed.point;
  State              state = seed.state;
  vnl_matrix<double> p(seed.covariance);
  double             fa = seed.fa;
  double             fa2 = seed.fa2;
  double             dNormMSE = 0; // no error at the seed
  double             trace = seed.trace;
  double             trace2 = seed.trace2;

  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

  vec_t m1, l1, m2, l2, m3, l3;
  m1 = seed.start_dir;

  // Tract the fiber.
  vnl_matrix<double> signal_tmp(_model->signal_dim(), 1);
  vnl_matrix<double> state_tmp(_model->state_dim(), 1);

  int stepnr = 0;

  while( true )
    {
    ++stepnr;

    Step3T(thread_id, x, m1, l1, m2, l2, m3, l3, fa, fa2, state, p, dNormMSE, trace, trace2);

    // Check if we should abort following this fiber. We abort if we reach the
    // CSF, if FA or GA get too small, if the curvature get's too high or if
    // the fiber gets too long.

    vnl_vector_ref<double> state_vnl(state.size(), &state.front() );

    bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1;

    state_tmp.set_column(0, state_vnl);
    _model->H(state_tmp, signal_tmp);

    double ga = s2ga(signal_tmp);
    bool   in_csf = ga < _ga_min || fa < _fa_min;
    bool   is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain || in_csf
        || static_cast<int>(fiber.position.size() ) > _max_length  // Stop if the fiber is too long
        || is_curving )
      {

      break;

      }

    Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

    // Record branch if necessary.
    if( is_branching )
      {
      bool is_one = l1._[0] > l1._[1] && l1._[0] > l1._[2];
      is_one = is_one && l2fa(l1._[0], l1._[1], l1._[2]) > _fa_min;
      if( is_one )
        {
        bool add_m2 = false;
        bool add_m3 = false;

        bool is_two = l2._[0] > l2._[1] && l2._[0] > l2._[2];
        bool is_three = l3._[0] > l3._[1] && l3._[0] > l3._[2];
        is_two = is_two && l2fa(l2._[0], l2._[1], l2._[2]) > _fa_min;
        is_three = is_three && l2fa(l3._[0], l3._[1], l3._[2]) > _fa_min;

        bool is_branch1 =
          dot(m1, m2) < _cos_theta_min && dot(m1, m2) > _cos_theta_max;
        bool is_branch2 =
          dot(m1, m3) < _cos_theta_min && dot(m1, m3) > _cos_theta_max;
        bool is_branch3 =
          dot(m2, m3) < _cos_theta_min;

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
              if( l2fa(l2._[0], l2._[1], l2._[2]) >
                  l2fa(l3._[0], l3._[1], l3._[2]) )
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

          local_seed.covariance.set_size(state_dim, state_dim);
          local_seed.covariance = p;

          SwapState3T(local_seed.state, local_seed.covariance, 2);
          local_seed.point = x;
          local_seed.start_dir = m2;
          local_seed.fa = l2fa(l2._[0], l2._[1], l2._[2]);
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
          local_seed.covariance.set_size(state_dim, state_dim);
          local_seed.covariance = p;
          SwapState3T(local_seed.state, local_seed.covariance, 3);
          local_seed.point = x;
          local_seed.start_dir = m3;
          local_seed.fa = l2fa(l3._[0], l3._[1], l3._[2]);
          }
        }
      }

    }
}

void Tractography::Follow2T(const int thread_id,
                            const size_t seed_index,
                            const SeedPointInfo& seed,
                            UKFFiber& fiber,
                            bool is_branching,
                            std::vector<SeedPointInfo>& branching_seeds,
                            std::vector<BranchingSeedAffiliation>& branching_seed_affiliation)
{

  // Unpack the seed information.
  vec_t x = seed.point;   // NOTICE that the x here is in ijk coordinate system
  State state = seed.state;

  vnl_matrix<double> p(seed.covariance);

  double fa = seed.fa;
  double fa2 = seed.fa2;
  double dNormMSE = 0; // no error at the seed
  double trace = seed.trace;
  double trace2 = seed.trace2;

  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2); // writes state to fiber.state

  vec_t m1, l1, m2, l2;
  m1 = seed.start_dir;

  // Track the fiber.
  vnl_matrix<double> signal_tmp(_model->signal_dim(), 1);
  vnl_matrix<double> state_tmp(_model->state_dim(), 1);

  int stepnr = 0;

  // useful for debugging
//   std::ofstream stateFile;
//   stateFile.open("states.txt", std::ios::app);

  while( true )
    {
    ++stepnr;

    Step2T(thread_id, x, m1, l1, m2, l2, fa, fa2, state, p, dNormMSE, trace, trace2);

    // Check if we should abort following this fiber. We abort if we reach the
    // CSF, if FA or GA get too small, if the curvature get's too high or if
    // the fiber gets too long.
    bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1; // is this 0.1 correct? yes

    // after here state doesnt change until next step.

    vnl_vector_ref<double> state_vnl(state.size(), &state.front() );

    // stateFile << state_vnl << std::endl;

    state_tmp.set_column(0, state_vnl);

    _model->H(state_tmp, signal_tmp); // signal_tmp is written, but only used to calculate ga

    double ga = s2ga(signal_tmp);
    bool   in_csf = ga < _ga_min || fa < _fa_min;
    // bool in_csf = ga < _ga_min || fa < _fa_min ;
    bool is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain
        || in_csf
        || static_cast<int>(fiber.position.size() ) > _max_length  // Stop when the fiber is too long
        || is_curving )
      {
      break;
      }

    Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

    // Record branch if necessary.
    if( is_branching )
      {
      bool is_two = l1._[0] > l1._[1] && l1._[0] > l1._[2] &&
        l2._[0] > l2._[1] && l2._[0] > l2._[2];
      is_two = is_two && l2fa(l1._[0], l1._[1], l1._[2]) > _fa_min &&
        l2fa(l2._[0], l2._[1], l2._[2]) > _fa_min;
      double theta = dot(m1, m2);
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
        local_seed.covariance.set_size(state_dim, state_dim);
        local_seed.covariance = p;
        SwapState2T(local_seed.state, local_seed.covariance);
        local_seed.point = x;
        local_seed.start_dir = m2;
        local_seed.fa = l2fa(l2._[0], l2._[1], l2._[2]);
        }
      }
    }

//   stateFile.close();
}

// Also read the comments to Follow2T above, it's documented better than this
// function here.
void Tractography::Follow1T(const int thread_id,
                            const SeedPointInfo& seed,
                            UKFFiber& fiber)
{

  assert(_model->signal_dim() == _signal_data->GetSignalDimension() * 2);

  vec_t x = seed.point;
  State state = seed.state;

  // DEBUG
//   std::cout << "seed state:\n";
//   for (int i=0;i<state.size();++i) {
//     std::cout << state[i] << " ";
//   }
//   std::cout << std::endl;

  vnl_matrix<double> p(seed.covariance);

  double fa = seed.fa;
  double fa2 = seed.fa2; // just needed for record
  double trace = seed.trace;
  double trace2 = seed.trace2;

  double dNormMSE = 0; // no error at the seed
  // Record start point.
  Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

  // Tract the fiber.
  vnl_matrix<double> signal_tmp(_model->signal_dim(), 1);
  vnl_matrix<double> state_tmp(_model->state_dim(), 1);

  int stepnr = 0;

  while( true )
    {
    ++stepnr;

    Step1T(thread_id, x, fa, state, p, dNormMSE, trace);

    // Terminate if off brain or in CSF.
    bool is_brain = _signal_data->Interp3ScalarMask(x) > 0.1; // x is the seed point

    vnl_vector_ref<double> state_vnl(state.size(), &state.front() );

    state_tmp.set_column(0, state_vnl);

    _model->H(state_tmp, signal_tmp);

    double ga = s2ga(signal_tmp);
    bool   in_csf = ga < _ga_min || fa < _fa_min;
    bool   is_curving = curve_radius(fiber.position) < _min_radius;

    if( !is_brain
        || in_csf
        || static_cast<int>(fiber.position.size() ) > _max_length  // Stop when fiber is too long
        || is_curving )
      {

      break;

      }

    Record(x, fa, fa2, state, p, fiber, dNormMSE, trace, trace2);

    }

}

void Tractography::Step3T(const int thread_id,
                          vec_t& x,
                          vec_t& m1,
                          vec_t& l1,
                          vec_t& m2,
                          vec_t& l2,
                          vec_t& m3,
                          vec_t& l3,
                          double& fa,
                          double& fa2,
                          State& state,
                          vnl_matrix<double>& covariance,
                          double& dNormMSE,
                          double& trace,
                          double& trace2
                          )
{

  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
         static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );
  State state_new(_model->state_dim() );

  vnl_matrix<double> covariance_new(_model->state_dim(), _model->state_dim() );

  // Use the Unscented Kalman Filter to get the next estimate.
  std::vector<double> signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);
  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  vec_t old_dir = m1;

  _model->State2Tensor3T(state, old_dir, m1, l1, m2, l2, m3, l3);
  trace = l1._[0] + l1._[1] + l1._[2];
  trace2 = l2._[0] + l2._[1] + l2._[2];

  double dot1 = dot(m1, old_dir);
  double dot2 = dot(m2, old_dir);
  double dot3 = dot(m3, old_dir);
  if( dot1 < dot2 && dot3 < dot2 )
    {
    // Switch dirs and lambdas.
    vec_t tmp = m1;
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
    vec_t tmp = m1;
    m1 = m3;
    m3 = tmp;
    tmp = l1;
    l1 = l3;
    l3 = tmp;

    // Swap state.
    SwapState3T(state, covariance, 3);
    }

  // Update FA. If the first lamba is not the largest anymore the FA is set to
  // 0 what will lead to abortion in the tractography loop.
  if( l1._[0] < l1._[1] || l1._[0] < l1._[2] )
    {
    fa = 0.0;
    }
  else
    {
    fa = l2fa(l1._[0], l1._[1], l1._[2]);
    fa2 = l2fa(l2._[0], l2._[1], l2._[2]);
    }

  vec_t voxel = _signal_data->voxel();

  // CB: Bug corrected, dir._[i] should be divided by voxel._[i]
  vec_t dx = make_vec(m1._[2] / voxel._[0],
                      m1._[1] / voxel._[1],
                      m1._[0] / voxel._[2]);
  x = x + dx * _stepLength;

}

void Tractography::Step2T(const int thread_id,
                          vec_t& x,
                          vec_t& m1,
                          vec_t& l1,
                          vec_t& m2,
                          vec_t& l2,
                          double& fa,
                          double& fa2,
                          State& state,
                          vnl_matrix<double>& covariance,
                          double& dNormMSE,
                          double& trace,
                          double& trace2
                          )
{

  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
         static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );

  State              state_new(_model->state_dim() );
  vnl_matrix<double> covariance_new(_model->state_dim(), _model->state_dim() );
  covariance_new.fill(0);

  // Use the Unscented Kalman Filter to get the next estimate.
  std::vector<double> signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);

  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  const vec_t old_dir = m1;   // Direction in last step

  _model->State2Tensor2T(state, old_dir, m1, l1, m2, l2);   // The returned m1 and m2 are unit vector here
  trace = l1._[0] + l1._[1] + l1._[2];
  trace2 = l2._[0] + l2._[1] + l2._[2];

  const double fa_tensor_1 = l2fa(l1._[0], l1._[1], l1._[2]);
  const double fa_tensor_2 = l2fa(l2._[0], l2._[1], l2._[2]);

  const double tensor_angle = std::acos(dot(m1, m2) ) * (180 / M_PI);

  if( dot(m1, old_dir) < dot(m2, old_dir) )
    {
    // Switch dirs and lambdas.
    vec_t tmp = m1;
    m1 = m2;
    m2 = tmp;
    tmp = l1;
    l1 = l2;
    l2 = tmp;

    vnl_matrix<double> old = covariance;

    SwapState2T(state, covariance);   // Swap the two tensors

    // DEBUG:
    // std::cout << "iHappened: " << std::acos(dot(m1,m2) )*(180/M_PI) << "\n";

    }

  if( tensor_angle <= 20 && std::min(fa_tensor_1, fa_tensor_2) <= 0.2 )
    {
    if( fa_tensor_1 > 0.2 )
      {
      // do nothing
      // i.e. keep m1 as principal direction
      }
    else
      {
      // switch directions, note: FA will be calculated alter
      vec_t tmp = m1;
      m1 = m2;
      m2 = tmp;
      tmp = l1;
      l1 = l2;
      l2 = tmp;

      vnl_matrix<double> old = covariance;

      SwapState2T(state, covariance);   // Swap the two tensors

      // DEBUG:
      // std::cout << "iHappened: " << std::acos(dot(m1,m2) )*(180/M_PI) << "\n";
      }
    }

  // Update FA. If the first lamba is not the largest anymore the FA is set to
  // 0 what will lead to abortion in the tractography loop.
  if( l1._[0] < l1._[1] || l1._[0] < l1._[2] )
    {
    fa = 0.0;
    fa2 = 0.0;
    }
  else
    {
    fa = l2fa(l1._[0], l1._[1], l1._[2]);
    fa2 = l2fa(l2._[0], l2._[1], l2._[2]);
    }

  vec_t dir = m1;     // The dir is a unit vector in ijk coordinate system indicating the direction of step

  vec_t voxel = _signal_data->voxel();

  vec_t dx = make_vec(dir._[2] / voxel._[0],  // By dividing by the voxel size, it's guaranteed that the step
                                              // represented by dx is 1mm in RAS coordinate system, no matter whether
                                              // the voxel is isotropic or not
                      dir._[1] / voxel._[1],  // The value is scaled back during the ijk->RAS transformation when
                                              // outputted
                      dir._[0] / voxel._[2]);

  x = x + dx * _stepLength; // The x here is in ijk coordinate system.
  // NOTICE that the coordinate order of x is in reverse order with respect to the axis order in the original signal
  // file.
  // This coordinate order is filpped back during output
  // The step length is in World space
  // exit(1);

}

void Tractography::Step1T(const int thread_id,
                          vec_t& x,
                          double& fa,
                          State& state,
                          vnl_matrix<double>& covariance,
                          double& dNormMSE,
                          double& trace
                          )
{

  assert(static_cast<int>(covariance.cols() ) == _model->state_dim() &&
         static_cast<int>(covariance.rows() ) == _model->state_dim() );
  assert(static_cast<int>(state.size() ) == _model->state_dim() );
  State              state_new(_model->state_dim() );
  vnl_matrix<double> covariance_new(_model->state_dim(), _model->state_dim() );

  std::vector<double> signal(_signal_data->GetSignalDimension() * 2);
  _signal_data->Interp3Signal(x, signal);

  _ukf[thread_id]->Filter(state, covariance, signal, state_new, covariance_new, dNormMSE);

  state = state_new;
  covariance = covariance_new;

  vec_t dir, l;
  _model->State2Tensor1T(state, dir, l);

  trace = l._[0] + l._[1] + l._[2];

  // Update FA. If the first lamba is not the largest anymore the FA is set to
  // 0 what will lead to abortion in the tractography loop.
  if( l._[0] < l._[1] || l._[0] < l._[2] )
    {
    fa = 0.0;
    }
  else
    {
    fa = l2fa(l._[0], l._[1], l._[2]);
    }

  vec_t voxel = _signal_data->voxel();

  vec_t dx = make_vec(dir._[2] / voxel._[0],
                      dir._[1] / voxel._[1],
                      dir._[0] / voxel._[2]);
  x = x + dx * _stepLength;

}

void Tractography::SwapState3T(State& state,
                               vnl_matrix<double>& covariance,
                               int i)
{

  // This function is only for 3T.
  assert(i == 2 || i == 3);

  vnl_vector_ref<double> state_vnl(state.size(), &state.front() );

  int                state_dim = _model->state_dim();
  vnl_matrix<double> tmp(state_dim, state_dim);
  state_dim /= 3;
  assert(state_dim == 5 || state_dim == 6);
  --i;
  int j = i == 1 ? 2 : 1;
  i *= state_dim;
  j *= state_dim;

  tmp.fill(0);
  tmp = covariance;
  covariance.update(tmp.extract(state_dim, state_dim, 0, 0), i, i);
  covariance.update(tmp.extract(state_dim, state_dim, i, i), 0, 0);
  covariance.update(tmp.extract(state_dim, state_dim, i, 0), 0, i);
  covariance.update(tmp.extract(state_dim, state_dim, 0, i), i, 0);

  covariance.update(tmp.extract(state_dim, state_dim, j, 0), j, i);
  covariance.update(tmp.extract(state_dim, state_dim, j, i), j, 0);
  covariance.update(tmp.extract(state_dim, state_dim, 0, j), i, j);
  covariance.update(tmp.extract(state_dim, state_dim, i, j), 0, j);

  // Swap the state.

  vnl_vector<double> tmp_vec(state_dim * 3);
  tmp_vec = state_vnl;

  state_vnl.update(tmp_vec.extract(state_dim, 0), i);
  state_vnl.update(tmp_vec.extract(state_dim, i), 0);

}

void Tractography::SwapState2T( State& state,
                                vnl_matrix<double>& covariance)
{
  // This function is only for 2T.

  vnl_vector_ref<double> state_vnl(state.size(), &state.front() );

  int state_dim = _model->state_dim();

  vnl_matrix<double> tmp(state_dim, state_dim);
  bool               bUnevenState = false;

  if( state_dim % 2 != 0 )
    {
    bUnevenState = true;                       // there is a weight term in the end of the state
    }
  state_dim = state_dim >> 1; // for uneven state (fw) rounds down, thats good

  tmp.fill(0);
  tmp = covariance;

  covariance.update(tmp.extract(state_dim, state_dim, 0, 0), state_dim, state_dim);
  covariance.update(tmp.extract(state_dim, state_dim, state_dim, state_dim), 0, 0);
  covariance.update(tmp.extract(state_dim, state_dim, state_dim, 0), 0, state_dim);
  covariance.update(tmp.extract(state_dim, state_dim, 0, state_dim), state_dim, 0);

  if( bUnevenState )   // change covariances of weights and state so they match the state again
    {
    covariance.update(tmp.extract(1, state_dim, state_dim * 2, 0), state_dim * 2, state_dim);
    covariance.update(tmp.extract(1, state_dim, state_dim * 2, state_dim), state_dim * 2, 0);

    covariance.update(tmp.extract(state_dim, 1, 0, state_dim * 2), state_dim, state_dim * 2);
    covariance.update(tmp.extract(state_dim, 1, state_dim, state_dim * 2), 0, state_dim * 2);
    }

  // Swap the state.

  vnl_vector<double> tmp_vec(state_dim * 2);
  tmp_vec = state_vnl;
  state_vnl.update(tmp_vec.extract(state_dim, 0), state_dim);
  state_vnl.update(tmp_vec.extract(state_dim, state_dim), 0);

}

void Tractography::Record(const vec_t& x, double fa, double fa2, const State& state,
                          const vnl_matrix<double> p,
                          UKFFiber& fiber, double dNormMSE, double trace, double trace2)
{

  assert(_model->state_dim() == static_cast<int>(state.size() ) );
  assert(p.rows() == static_cast<unsigned int>(state.size() ) &&
         p.cols() == static_cast<unsigned int>(state.size() ) );

  // std::cout << "x: " << x._[0] << " " << x._[1] << " " << x._[2] << std::endl;

  fiber.position.push_back(x);
  fiber.norm.push_back(p.array_two_norm() );

  if( _record_nmse )
    {
    fiber.normMSE.push_back(dNormMSE);
    }

  if( _record_trace )
    {
    fiber.trace.push_back(trace);
    if( _num_tensors >= 2 )
      {
      fiber.trace2.push_back(trace2);
      }
    }

  if( _record_fa )
    {
    fiber.fa.push_back(fa);
    if( _num_tensors >= 2 )
      {
      fiber.fa2.push_back(fa2);
      }
    }

  if( _record_free_water )
    {
    double fw = 1 - state[_nPosFreeWater];
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
  if( state.size() == 5 || state.size() == 10 || state.size() == 15 )   // i.e. simple model
    { // Normalize direction before storing it;
    State store_state(state);
    vec_t dir = make_vec(store_state[0], store_state[1], store_state[2]);
    dir /= norm(dir);
    store_state[0] = dir._[0];
    store_state[1] = dir._[1];
    store_state[2] = dir._[2];

    if( state.size() == 10 )
      {
      dir = make_vec(store_state[5], store_state[6], store_state[7]);
      dir /= norm(dir);
      store_state[5] = dir._[0];
      store_state[6] = dir._[1];
      store_state[7] = dir._[2];
      }
    if( state.size() == 15 )
      {
      dir = make_vec(store_state[10], store_state[11], store_state[12]);
      dir /= norm(dir);
      store_state[10] = dir._[0];
      store_state[11] = dir._[1];
      store_state[12] = dir._[2];
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
