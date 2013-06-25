/**
 * \file vtk_writer.h
 * \brief Fiber writing functionality
 *
 * Contains Class definition of the VtkWriter, used for storing the fiber vector to a .vtk file
 * \todo Could be done more elegantly with VTK
*/

#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include <string>
#include <vector>
#include "linalg.h"
#include "tractography.h"

struct Fiber;

class ISignalData;

/**
 * \class VtkWriter
 * \brief Class that allows to write a bunch of fibers to a .vtk file
*/
class VtkWriter
{
public:
  /**
   * \brief Constructor
  */
  VtkWriter(const ISignalData *signal_data, Tractography::model_type filter_model_type,  bool write_tensors);

  /** Destructor */
  virtual ~VtkWriter()
  {
  }

  /**
   * \brief Writes the fibers to the VTK file and attaches the selected values to the fiber
   * \param[in] file_name The path of the output fiber file (*.vtk)
   * \param[in] tractsWithSecondTensor File path for the fibers generated with the second tensor
   *                                   This one is optional.
   * \param[in] store_glyphs Write glyphs (i.e. main tensor directions) to a file named glyphs_{tracts}.
  */
  bool Write(const std::string& file_name, const std::string & tractsWithSecondTensor, const std::vector<Fiber>& fibers,
             bool write_state, bool store_glyphs);

  /** Write the glyphs (i.e. main tensor directions) to  a file named glyphs_{tracts}. */
  bool WriteGlyphs(const std::string& file_name, const std::vector<Fiber>& fibers);

  /** Sets the variable that toggles the transform from ijk to RAS before writing the fiber to VTK. */
  void set_transform_position(bool transform_position)
  {
    _transform_position = transform_position;
  }

protected:
  /**
   * Writes a point to the output stream, and performs the ijk-RAS transform if set to do so.
   * Also this function makes sure that the output VTK file is well formatted.
  */
  void WritePoint(const vec_t& point, std::ofstream& output, int& counter);

  /**
   * Writes the fibers and all values attached to them to a VTK file
  */
  void writeFibersAndTensors(std::ofstream & output, const std::vector<Fiber>& fibers, const int tensorNumber);

  /**
   * \brief Reconstructs the tensor from the state for each case
   * \param[out] D The calculated diffusion tensor
   * \todo I think there is something wrong with choosing a orthonormal basis for the tensor
  */
  void State2Tensor(State & state, mat_t & D, int tensorNumber);

  /** The diffusion weighted signal data */
  const ISignalData *_signal_data;

  /** If set to true, WritePoint will convert the points from ijk to RAS space before writing */
  bool _transform_position;

  /** What model was the tractography performed with */
  const Tractography::model_type _filter_model_type;

  /** Scaling of the glyphs */
  const double _scale_glyphs;

  /** Whether to attach the tensors to the fiber */
  bool _write_tensors;

  // Positions of the eigenvalues, directions, or angles, respectively, in the state.
  int _p_l1, _p_l2, _p_l3;
  int _p_m1, _p_m2, _p_m3;
  int _p_psi, _p_theta, _p_phi;

  /** How many indeces in the state are used up by 1 tensor */
  int _tensor_space;

  /** The number of tensors used */
  int _num_tensors;

  /** Was the full or simplified tensor model used? */
  bool _full;

  /** Additional scaling of the eigenvalues before writing */
  const double _eigenScaleFactor;

  /** Transformation matrix from ijk-RAS with voxel size normalized out */
  mat_t _sizeFreeI2R;
};

#endif  // VTK_WRITER_H_
