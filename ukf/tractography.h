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
#include "ukf_types.h"

class ISignalData;

struct UKFSettings {
  bool record_fa;
  bool record_nmse;
  bool record_trace;
  bool record_state;
  bool record_cov;
  bool record_free_water;
  bool record_tensors;
  bool record_Vic;
  bool record_kappa;
  bool record_Viso;
  bool transform_position;
  bool store_glyphs;
  bool branches_only;
  ukfPrecisionType fa_min;
  ukfPrecisionType mean_signal_min;
  ukfPrecisionType seeding_threshold;
  int num_tensors;
  int seeds_per_voxel;
  ukfPrecisionType min_branching_angle;
  ukfPrecisionType max_branching_angle;
  bool is_full_model;
  bool free_water;
  bool noddi;
  ukfPrecisionType stepLength;
  ukfPrecisionType recordLength;
  ukfPrecisionType maxHalfFiberLength;
  std::vector<int> labels;

  ukfPrecisionType p0;
  ukfPrecisionType sigma_signal;
  ukfPrecisionType sigma_mask;
  ukfPrecisionType min_radius;
  ukfPrecisionType full_brain_mean_signal_min;
  size_t num_threads;
};

/**
 * \class Tractography
 * \brief This class performs the tractography and saves each step.
*/
class Tractography
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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
  Tractography(UKFSettings settings,
               FilterModel *model,
               model_type filter_model_type,
               const std::string output_file,
               const std::string output_file_with_second_tensor);

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

  /** \breif Performs the tractography
      \return true if files written successfully, else false
  */
  bool Run();

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
  void SetWriteCompressed(bool wb) { this->_writeCompressed = wb; }
private:
  /**
   * Calculate six tensor coefficients by solving B * d = log(s), where d are
   * tensor coefficients, B is gradient weighting, s is signal.
  */
  void UnpackTensor(const ukfVectorType& b, const stdVec_t& u, stdEigVec_t& s,
                    stdEigVec_t& ret);

  /**
  * Creates necessary variable for noddi
  */
  void createProtocol(const ukfVectorType& b, ukfVectorType& gradientStrength,
                                  ukfVectorType& pulseSeparation);

  /** One step along the fiber for the 3-tensor case. */
  void Step3T(const int thread_id, vec3_t& x, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2, vec3_t& m3, vec3_t& l3,
              ukfPrecisionType& fa, ukfPrecisionType& fa2, State& state, ukfMatrixType& covariance, ukfPrecisionType& dNormMSE, ukfPrecisionType& trace,
              ukfPrecisionType& trace2);

  /** One step along the fiber for the 2-tensor case. */
  void Step2T(const int thread_id, vec3_t& x, vec3_t& m1, vec3_t& l1, vec3_t& m2, vec3_t& l2, ukfPrecisionType& fa, ukfPrecisionType& fa2,
              State& state, ukfMatrixType& covariance, ukfPrecisionType& dNormMSE, ukfPrecisionType& trace, ukfPrecisionType& trace2);

  /** One step along the fiber for the 1-tensor case. */
  void Step1T(const int thread_id, vec3_t& x, ukfPrecisionType& fa, State& state, ukfMatrixType& covariance, ukfPrecisionType& dNormMSE,
              ukfPrecisionType& trace);

  /**
   * Swaps the first tensor with the i-th tensor in state and covariance matrix for the 3 Tensor case.
   * This is used when the main direction of the tractography 'switches' tensor.
  */
  void SwapState3T(State& state, ukfMatrixType& covariance, int i);

  /**
   * Swap the tensors in the state and covariance matrix for the 2-tensor case. This is used when the
   * principal direction of the minor tensor has more weight than the one of the current tensor.
  */
  void SwapState2T(State& state, ukfMatrixType& covariance);

  /**
   * Saves one point along the fiber so that everything can be written to a
   * file at the end.
  */
  void Record(const vec3_t& x, const ukfPrecisionType fa, const ukfPrecisionType fa2,
              const State& state, const ukfMatrixType p, UKFFiber& fiber,
              const ukfPrecisionType dNormMSE, const ukfPrecisionType trace, const ukfPrecisionType trace2);

  /**  Reserving fiber array memory so as to avoid resizing at every step*/
  void FiberReserve(UKFFiber& fiber, int fiber_size);

  /** Vector of Pointers to Unscented Kalaman Filters. One for each thread. */
  std::vector<UnscentedKalmanFilter *> _ukf;

  /** Pointer to generic filter model */
  FilterModel *_model;

  /** Type of the filter model */
  model_type _filter_model_type;

  /** Output file for tracts generated with first tensor */
  const std::string _output_file;
  /** Output file for tracts generated with second tensor */
  const std::string _output_file_with_second_tensor;

  /** Pointer to generic diffusion data */
  ISignalData *_signal_data;

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
  // Noddi Model parameters
  /** Switch for attaching the Vic to the fiber at each point of the tractography */
  const bool _record_Vic;
  /** Switch for attaching the kappa to the fiber at each point of the tractography */
  const bool _record_kappa;
  /** Switch for attaching the Viso to the fiber at each point of the tractography */
  const bool _record_Viso;
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
  const ukfPrecisionType _p0;
  const ukfPrecisionType _sigma_signal;
  const ukfPrecisionType _sigma_mask;
  const ukfPrecisionType _min_radius;

  /** Maximal number of points in the tract */
  const int _max_length;
  bool      _full_brain;
  bool _noddi;
  ukfVectorType _gradientStrength, _pulseSeparation;
  /** Index of the weight in the state for the free water cases */
  int _nPosFreeWater;

  // Parameters for the tractography
  const ukfPrecisionType           _fa_min;
  const ukfPrecisionType           _mean_signal_min;
  const ukfPrecisionType           _seeding_threshold;
  const int              _num_tensors;
  const int              _seeds_per_voxel;
  ukfPrecisionType                 _cos_theta_min;
  ukfPrecisionType                 _cos_theta_max;
  const bool             _is_full_model;
  const bool             _free_water;
  const ukfPrecisionType           _stepLength;
  const int                 _steps_per_record;
  const std::vector<int> _labels;

  bool _writeBinary;
  bool _writeCompressed;
  // Threading control
  const int _num_threads;
};

#endif  // TRACTOGRAPHY_H_
