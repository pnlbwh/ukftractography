/**
 * \file tractography.h
 * \brief Contains the Class Tractography, which contains the functions that deal with the
 * actual tracing of the fibers for each model
*/
#ifndef TRACTOGRAPHY_H_
#define TRACTOGRAPHY_H_

#include <string>
#include <vector>
#include "ukffiber.h"
#include "seed.h"

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_vector_ref.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_qr.h>

class ISignalData;

/**
 * \class Tractography
 * \brief This class performs the tractography and saves each step.
*/
class Tractography
{

public:

  /** Defines a type to specify the model */
  enum model_type { _1T,
                    _1T_FW,
                    _1T_FULL,
                    _1T_FW_FULL,
                    _2T,
                    _2T_FW,
                    _2T_FULL,
                    _2T_FW_FULL,
                    _3T,
                    _3T_FULL };

  /** Constructor, is called from main.cc where all parameters are defined. */
  Tractography(FilterModel *model, model_type filter_model_type,

               const std::string& output_file, const std::string & output_file_with_second_tensor,
               const bool record_fa, const bool record_nmse, const bool record_trace,
               const bool record_state, const bool record_cov, const bool record_free_water, const bool record_tensors,
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
               );

  /** Destructor */
  ~Tractography();

  /**
   * Load the files that contain the DWI signal, the seeds and a mask
   * defining the volume of the brain.
  */
  bool LoadFiles(const std::string& data_file, const std::string& seed_file, const std::string& mask_file,
                 const bool normalized_DWI_data, const bool output_normalized_DWI_data);

  /**
   * Creates the seeds and initilizes them by finding the tensor directions,
   * eigenvalues and Euler Angles. This also sets the initial state and
   * covariance.
  */
  void Init(std::vector<SeedPointInfo>& seed_infos);

  /** Performs the tractography */
  void Run();

  /**
   * Follows one seed point for the 3 Tensor case
  */
  void Follow3T(const int thread_id, const size_t seed_index, const SeedPointInfo& seed, UKFFiber& fiber,
                bool is_branching, std::vector<SeedPointInfo>& branching_seeds,
                std::vector<BranchingSeedAffiliation>& branching_seed_affiliation);

  /**
   * Follows one seed point for the 2 Tensor case
  */
  void Follow2T(const int thread_id, const size_t seed_index, const SeedPointInfo& seed, UKFFiber& fiber,
                bool is_branching, std::vector<SeedPointInfo>& branching_seeds,
                std::vector<BranchingSeedAffiliation>& branching_seed_affiliation);

  /**
   * Follows one seed point for the 1 Tensor case
  */
  void Follow1T(const int thread_id, const SeedPointInfo& seed, UKFFiber& fiber);

  void SetWriteBinary(bool wb) { this->_writeBinary = wb; }
private:
  /**
   * Calculate six tensor coefficients by solving B * d = log(s), where d are
   * tensor coefficients, B is gradient weighting, s is signal.
  */
  void UnpackTensor(const std::vector<double>& b, const std::vector<vec_t>& u, std::vector<std::vector<double> >& s,
                    std::vector<std::vector<double> >& ret);

  /** One step along the fiber for the 3-tensor case. */
  void Step3T(const int thread_id, vec_t& x, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2, vec_t& m3, vec_t& l3,
              double& fa, double& fa2, State& state, vnl_matrix<double>& covariance, double& dNormMSE, double& trace,
              double& trace2);

  /** One step along the fiber for the 2-tensor case. */
  void Step2T(const int thread_id, vec_t& x, vec_t& m1, vec_t& l1, vec_t& m2, vec_t& l2, double& fa, double& fa2,
              State& state, vnl_matrix<double>& covariance, double& dNormMSE, double& trace, double& trace2);

  /** One step along the fiber for the 1-tensor case. */
  void Step1T(const int thread_id, vec_t& x, double& fa, State& state, vnl_matrix<double>& covariance, double& dNormMSE,
              double& trace);

  /**
   * Swaps the first tensor with the i-th tensor in state and covariance matrix for the 3 Tensor case.
   * This is used when the main direction of the tractography 'switches' tensor.
  */
  void SwapState3T(State& state, vnl_matrix<double>& covariance, int i);

  /**
   * Swap the tensors in the state and covariance matrix for the 2-tensor case. This is used when the
   * principal direction of the minor tensor has more weight than the one of the current tensor.
  */
  void SwapState2T(State& state, vnl_matrix<double>& covariance);

  /**
   * Saves one point along the fiber so that everything can be written to a
   * file at the end.
  */
  void Record(const vec_t& x, double fa, double fa2, const State& state, const vnl_matrix<double> p, UKFFiber& fiber,
              double dNormMSE, double trace, double trace2);

  /** Vector of Pointers to Unscented Kalaman Filters. One for each thread. */
  std::vector<UnscentedKalmanFilter *> _ukf;

  /** Pointer to generic diffusion data */
  ISignalData *_signal_data;
  /** Pointer to generic filter model */
  FilterModel *_model;

  /** Type of the filter model */
  model_type _filter_model_type;

  /** Output file for tracts generated with first tensor */
  const std::string _output_file;
  /** Output file for tracts generated with second tensor */
  const std::string _output_file_with_second_tensor;
  /** Switch for attaching the FA value to the fiber at each point of the tractography */
  const bool _record_fa;
  /**
   * Switch for attaching the normalized mean squared error of the reconstructed signal to the real signal
   * to the fiber at each point of the tractography
  */
  const bool _record_nmse;
  /** Switch for attaching the trace to the fiber at each point of the tractography */
  const bool _record_trace;
  /** Switch for attaching the state to the fiber at each point of the tractography */
  const bool _record_state;
  /** Switch for attaching the covariance to the fiber at each point of the tractography */
  const bool _record_cov;
  /** Switch for attaching the free water percentage to the fiber at each point of the tractography */
  const bool _record_free_water;
  /**
   * Switch for attaching the diffusion tensors to the fiber at each point of the tractography.
   * This is important for visualizing the fiber properties in Slicer.
  */
  const bool _record_tensors;
  /**
   * Wheather to transform the points back to RAS-space before writing the VTK or not.
  */
  const bool _transform_position;
  /** Attach the glyphs to the VTK file */
  const bool _store_glyphs;
  /** To output branches only */
  bool _branches_only;

  // Internal parameters
  bool         _is_branching;
  const double _p0;
  const double _sigma_signal;
  const double _sigma_mask;
  const double _min_radius;
  const double _full_brain_ga_min;
  /** Maximal number of points in the tract */
  const int _max_length;
  bool      _full_brain;
  /** Index of the weight in the state for the free water cases */
  int _nPosFreeWater;

  // Parameters for the tractography
  const double           _fa_min;
  const double           _ga_min;
  const double           _seedFALimit;
  const int              _num_tensors;
  const int              _seeds_per_voxel;
  double                 _cos_theta_min;
  double                 _cos_theta_max;
  const bool             _is_full_model;
  const bool             _free_water;
  const double           _stepLength;
  const std::vector<int> _labels;

  bool _writeBinary;

  // Threading control
  const int _num_threads;
};

#endif  // TRACTOGRAPHY_H_
