/**
 * \file vtk_writer.cc
 * \brief implementation of vtk_writer.h
 * \todo The case differentiation in the beginning is very hackish..
*/

#include "vtk_writer.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "ISignalData.h"
#include "utilities.h"
#include "ukffiber.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

/** Historical VTK File formats are VERY sensitive to their format **/
/** http://www.vtk.org/VTK/img/file-formats.pdf **/
static const char * const VTK_LEGACY_FORMAT_HEADER = "# vtk DataFile Version 3.0";

VtkWriter::VtkWriter(const ISignalData *signal_data, Tractography::model_type filter_model_type, bool write_tensors) :
  _signal_data(signal_data),
  _transform_position(true),
  _filter_model_type(filter_model_type),
  _scale_glyphs(0.01),
  _write_tensors(write_tensors),
  _eigenScaleFactor(1),
  _writeBinary(false)
{

  if( filter_model_type == Tractography::_1T || filter_model_type == Tractography::_1T_FW )
    {
    _full = false;
    _p_m1 = 0;
    _p_m2 = 1;
    _p_m3 = 2;
    _p_l1 = 3;
    _p_l2 = 4;
    _p_l3 = 4;
    _num_tensors = 1;
    _tensor_space = 5;
    }
  else if( filter_model_type == Tractography::_1T_FULL || filter_model_type == Tractography::_1T_FW_FULL )
    {
    _full = true;
    _p_l1 = 3,
      _p_l2 = 4;
    _p_l3 = 5;
    _num_tensors = 1;
    _tensor_space = 6;
    }
  else if( filter_model_type == Tractography::_2T || filter_model_type == Tractography::_2T_FW )
    {
    _full = false;
    _p_m1 = 0;
    _p_m2 = 1;
    _p_m3 = 2;
    _p_l1 = 3;
    _p_l2 = 4;
    _p_l3 = 4;
    _num_tensors = 2;
    _tensor_space = 5;
    }
  else if( filter_model_type == Tractography::_2T_FULL || filter_model_type == Tractography::_2T_FW_FULL )
    {
    _full = true;
    _p_l1 = 3, _p_l2 = 4;
    _p_l3 = 5;
    _num_tensors = 2;
    _tensor_space = 6;
    }
  else if( filter_model_type == Tractography::_3T )
    {
    _full = false;
    _p_m1 = 0;
    _p_m2 = 1;
    _p_m3 = 2;
    _p_l1 = 3;
    _p_l2 = 4;
    _p_l3 = 4;
    _num_tensors = 3;
    _tensor_space = 5;
    }
  else if( filter_model_type == Tractography::_3T_FULL )
    {
    _full = true;
    _p_l1 = 3, _p_l2 = 4;
    _p_l3 = 5;
    _num_tensors = 3;
    _tensor_space = 6;
    }

  // this also upon initialization of writer, its the same for all
  vnl_matrix<double> i2r = _signal_data->i2r();
  vec_t              voxel = _signal_data->voxel();
  // Factor out the effect of voxel size
  _sizeFreeI2R = make_mat(
    i2r(0, 0) / voxel._[2], i2r(0, 1) / voxel._[1], i2r(0, 2) / voxel._[0],
    i2r(1, 0) / voxel._[2], i2r(1, 1) / voxel._[1], i2r(1, 2) / voxel._[0],
    i2r(2, 0) / voxel._[2], i2r(2, 1) / voxel._[1], i2r(2, 2) / voxel._[0]
    );

}


void VtkWriter::writeFibersAndTensors(std::ofstream & output, const std::vector<UKFFiber>& fibers, const int tensorNumber)
{
  // Write file version and identifier
  output << VTK_LEGACY_FORMAT_HEADER << std::endl;

  // Write header
  output << "Tracts generated with tensor " << tensorNumber << std::endl;

  // File format
  if(this->_writeBinary)
    {
    output << "BINARY" << std::endl;
    }
  else
    {
    output << "ASCII" << std::endl;
    }

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  // Write the points.
  output << "DATASET POLYDATA" << std::endl;

  output << "POINTS " << num_points << " float" << std::endl;

  // The points are written so that 3 points fits in one line
  int counter = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j )
      {
      WritePoint(fibers[i].position[j], output, counter);
      }
    }
  output << std::endl;

  // Write the lines.
  counter = 0;
  output << std::endl << "LINES " << num_fibers << " "
         << num_fibers + num_points << std::endl;
  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    if(!this->_writeBinary)
      {
      output << fiber_size;
      for( int j = 0; j < fiber_size; ++j )
        {
        output << " " << counter++;
        }
      output << std::endl;
      }
    else
      {
      this->WriteX<int,int>(output,fiber_size);
      for(int j = 0; j < fiber_size; ++j)
        {
        this->WriteX<int,int>(output,counter);
        }
      }
    }
  output << std::endl;

  /////
  // Dataset attribute part starts
  /////

  output << "POINT_DATA " << num_points << std::endl;

  /////
  // Tensor output
  /////
  typedef std::vector<double> State;

  if( _write_tensors )
    {
    for( int local_tensorNumber = 1; local_tensorNumber <= _num_tensors; ++local_tensorNumber )
      {
      output << "TENSORS Tensor" << local_tensorNumber << " float" << std::endl;
      for( int i = 0; i < num_fibers; i++ )
        {
        const int fiber_size = static_cast<int>(fibers[i].position.size() );
        for( int j = 0; j < fiber_size; ++j )
          {
          const State & state = fibers[i].state[j];
          mat_t         D;
          State2Tensor(state, D, local_tensorNumber);
          if(!this->_writeBinary)
            {
            output << D._[0] * _eigenScaleFactor << " "
                   << D._[1] * _eigenScaleFactor << " "
                   << D._[2] * _eigenScaleFactor << std::endl
                   << D._[3] * _eigenScaleFactor << " "
                   << D._[4] * _eigenScaleFactor << " "
                   << D._[5] * _eigenScaleFactor << std::endl
                   << D._[6] * _eigenScaleFactor << " "
                   << D._[7] * _eigenScaleFactor << " "
                   << D._[8] * _eigenScaleFactor << std::endl;
            }
          else
            {
            for(unsigned k = 0; k < 9; ++k)
              {
              this->WriteX<double,float>(output,D._[k]);
              }
            }
          }
        }
      }
    }
  output << std::endl;
}

bool VtkWriter::Write(const std::string& file_name, const std::string & tractsWithSecondTensor,
                      const std::vector<UKFFiber>& fibers,
                      bool write_state, bool store_glyphs)
{
  if( fibers.size() == 0 )
    {
    std::cout << "No fiber exists." << std::endl;
    return true;
    }

  std::ios_base::openmode openMode = std::ios_base::out;
  if(this->_writeBinary)
    {
    openMode |= std::ios_base::binary;
    }
  // Open file for writing.
  std::ofstream output(file_name.c_str(), openMode);
  if( !output.is_open() )
    {
    std::cout << "Failed to open " << file_name << "." << std::endl;
    return true;
    }
  std::cout << "Writing to " << file_name << "." << std::endl;
  writeFibersAndTensors(output, fibers, 1);

  if( !tractsWithSecondTensor.empty() )
    {
    std::ofstream output_withSecondTensor(tractsWithSecondTensor.c_str() );
    if( !output_withSecondTensor.is_open() )
      {
      std::cout << "Failed to open " << tractsWithSecondTensor << "." << std::endl;
      return true;
      }
    std::cout << "Writing to " << tractsWithSecondTensor << "." << std::endl;
    writeFibersAndTensors(output_withSecondTensor, fibers, 2);
    }

  // Glyph output
  if( store_glyphs )
    {
    assert(fibers[0].state.size() );
    std::stringstream ss;
    ss << file_name.substr(0, file_name.find_last_of(".") ) << "_glyphs"
       << ".vtk";
    if( WriteGlyphs(ss.str(), fibers) )
      {
      return true;
      }
    }

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  /////
  // Here starts a field data that may contain norm, FA, state and covariance
  /////
  int  fields = 1; // 1 for the norm.
  bool write_fa = fibers[0].fa.size();
  bool write_fa2 = fibers[0].fa2.size();
  bool write_cov = fibers[0].covariance.size();
  bool write_fw = fibers[0].free_water.size();
  bool write_nmse = fibers[0].normMSE.size();
  bool write_trace = fibers[0].trace.size();
  bool write_trace2 = fibers[0].trace2.size();

  if( write_fa )
    {
    ++fields;
    }
  if( write_fa2 )
    {
    ++fields;
    }
  if( write_state )
    {
    ++fields;
    }
  if( write_cov )
    {
    ++fields;
    }
  if( write_fw )
    {
    ++fields;
    }
  if( write_nmse )
    {
    ++fields;
    }
  if( write_trace )
    {
    ++fields;
    }
  if( write_trace2 )
    {
    ++fields;
    }
  output << "FIELD FieldData " << fields << std::endl;

  // TODO: write a function for writing the scalar fields

  // Write norm.
  int counter = 0;
  output << "norm 1 " << num_points << " float" << std::endl;
  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j )
      {
      if(!this->_writeBinary)
        {
        output << fibers[i].norm[j];
        ++counter;
        if( !(counter % 9) )
          {
          output << std::endl;
          }
        else if( counter < num_points )
          {
          output << " ";
          }
        }
      else
        {
        this->WriteX<double,float>(output,fibers[i].norm[j]);
        }
      }
    }
  if(this->_writeBinary || counter % 9 )
    {
    output << std::endl;
    }

  // Fractional anisotropy
  if( write_fa )
    {
    output << "FA 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].fa[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output, fibers[i].fa[j]);
          }
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  if( write_fa2 )
    {
    output << "FA2 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].fa2[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output, fibers[i].fa2[j]);
          }
        }
      }
    if(!this->_writeBinary &&  counter % 9 )
      {
      output << std::endl;
      }
    }

  // trace
  if( write_trace )
    {
    output << "Trace 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].trace[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output,fibers[i].trace[j]);
          }
        }
      }
    if(!this->_writeBinary &&  counter % 9 )
      {
      output << std::endl;
      }
    }

  // trace
  if( write_trace2 )
    {
    output << "Trace2 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].trace2[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output,fibers[i].trace2[j]);
          }
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  // Free water
  if( write_fw )
    {
    output << "FreeWater 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].free_water[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output,fibers[i].free_water[j]);
          }
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  // Normalized Mean Squared Error of the estimated signal to the real signal
  double nmse_sum = 0.0;
  if( write_nmse )
    {
    output << "NMSE 1 " << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        if(!this->_writeBinary)
          {
          output << fibers[i].normMSE[j];
          ++counter;
          if( !(counter % 9) )
            {
            output << std::endl;
            }
          else if( counter < num_points )
            {
            output << " ";
            }
          }
        else
          {
          this->WriteX<double,float>(output,fibers[i].normMSE[j]);
          }
        nmse_sum += fibers[i].normMSE[j];
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  std::cout << "nmse_avg=" << nmse_sum / counter << std::endl;

  // State output
  if( write_state )
    {
    int state_dim = fibers[0].state[0].size();
    output << "state " << state_dim << " " << num_points << " float"
           << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        for( int k = 0; k < state_dim; ++k )
          {
          if(!this->_writeBinary)
            {
            output << fibers[i].state[j][k];
            ++counter;
            if( !(counter % 9) )
              {
              output << std::endl;
              }
            else if( counter < num_points * state_dim )
              {
              output << " ";
              }
            }
          else
            {
            this->WriteX<double,float>(output,fibers[i].state[j][k]);
            }
          }
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  // Covariance info
  if( write_cov )
    {
    int state_dim = fibers[0].state[0].size();
    output << "covariance " << state_dim * (state_dim + 1) / 2 << " "
           << num_points << " float" << std::endl;
    counter = 0;
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        for( int a = 0; a < state_dim; ++a )
          {
          for( int b = a; b < state_dim; ++b )
            {
            if(!this->_writeBinary)
              {
              output << fibers[i].covariance[j] (a, b);
              ++counter;
              if( !(counter % 9) )
                {
                output << std::endl;
                }
              else if( counter < num_points * state_dim * (state_dim + 1) / 2 )
                {
                output << " ";
                }
              }
            else
              {
              this->WriteX<double,float>(output,fibers[i].covariance[j] (a, b));
              }
            }
          }
        }
      }
    if(this->_writeBinary || counter % 9 )
      {
      output << std::endl;
      }
    }

  output.close();

  return false;
}

bool VtkWriter::WriteGlyphs(const std::string& file_name,
                            const std::vector<UKFFiber>& fibers)
{
  if( fibers.size() == 0 )
    {
    std::cout << "No fibers existing." << std::endl;
    return true;
    }

  // Open file for writing.
  std::ofstream output(file_name.c_str() );
  if( !output.is_open() )
    {
    std::cout << "Failed to open " << file_name << "." << std::endl;
    return true;
    }

  std::cout << "Writing glyphs to " << file_name << "." << std::endl;

  // Write header information
  output << VTK_LEGACY_FORMAT_HEADER << std::endl;
  output << "TractographyGlyphs" << std::endl;
  output << "ASCII" << std::endl;

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  int num_tensors = fibers[0].state[0].size() / 5;

  const double scale = 0.5 * _scale_glyphs;

  // Write the points.
  output << "DATASET POLYDATA" << std::endl;
  output << "POINTS " << num_points * 2 * num_tensors << " float";
  int counter = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j )
      {
      vec_t        point = fibers[i].position[j];
      const State& state = fibers[i].state[j];

      // Get the directions.
      vec_t m1, m2, m3;
      if( state.size() == 5 )
        {
        m1 = make_vec(state[2], state[1], state[0]);
        m1 = state[3] / 100.0 * m1;
        }
      else if( state.size() == 6 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        double tmp = m1._[0];
        m1._[0] = m1._[2];
        m1._[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        }
      else if( state.size() == 10 )
        {
        m1 = make_vec(state[2], state[1], state[0]);
        m2 = make_vec(state[7], state[6], state[5]);
        m1 = state[3] / 100.0 * m1;
        m2 = state[8] / 100.0 * m2;
        }
      else if( state.size() == 12 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        m2 = rotation_main_dir(state[6], state[7], state[8]);
        double tmp = m1._[0];
        m1._[0] = m1._[2];
        m1._[2] = tmp;
        tmp = m2._[0];
        m2._[0] = m2._[2];
        m2._[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        m2 = state[9] / 100.0 * m2;
        }
      else if( state.size() == 15 )
        {
        m1 = make_vec(state[2], state[1], state[0]);
        m2 = make_vec(state[7], state[6], state[5]);
        m3 = make_vec(state[12], state[11], state[10]);
        m1 = state[3] / 100.0 * m1;
        m2 = state[8] / 100.0 * m2;
        m3 = state[13] / 100.0 * m3;
        }
      else if( state.size() == 18 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        m2 = rotation_main_dir(state[6], state[7], state[8]);
        m3 = rotation_main_dir(state[12], state[13], state[14]);
        double tmp = m1._[0];
        m1._[0] = m1._[2];
        m1._[2] = tmp;
        tmp = m2._[0];
        m2._[0] = m2._[2];
        m2._[2] = tmp;
        tmp = m3._[0];
        m3._[0] = m3._[2];
        m3._[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        m2 = state[9] / 100.0 * m2;
        m3 = state[15] / 100.0 * m3;
        }

      // Calculate the points. The glyphs are represented as two-point lines.
      vec_t pos1, pos2;
      if( num_tensors == 1 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);
        }
      else if( num_tensors == 2 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);
        }
      else if( num_tensors == 3 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);

        pos1 = point - scale * m3;
        pos2 = point + scale * m3;
        WritePoint(pos1, output, counter);
        WritePoint(pos2, output, counter);
        }
      }
    }
  if(!this->_writeBinary)
    {
    output << std::endl;
    }

  // Write the lines.
  output << std::endl << "LINES " << num_points * num_tensors << " "
         << num_points * num_tensors * 3 << std::endl;
  for( int i = 0; i < num_points * num_tensors; ++i )
    {
    output << "2 " << i * 2 << " " << i * 2 + 1 << std::endl;
    }

  output.close();

  return false;
}

void VtkWriter::WritePoint(const vec_t& point, std::ofstream& output,
                           int& counter)
{
  if(!this->_writeBinary)
    {
    // Always write nine floats onto one line.
    if( !(counter % 3) )
      {
      output << std::endl;
      }
    else
      {
      output << " ";
      }
    }
  // Transform the point into the correct coordinate system.
  vnl_vector<double> p(4);
  p[0] = point._[2];    // NOTICE the change of order here. Flips back to the original axis order
  p[1] = point._[1];
  p[2] = point._[0];

  if( _transform_position )
    {
    p[3] = 1.0;
    vnl_vector<double> p_new(4);
    p_new = _signal_data->i2r() * p;    // ijk->RAS transform
    // reverse order again, so that only one output statement is needed.
    p[2] = p_new[0];
    p[1] = p_new[1];
    p[0] = p_new[2];
    // output << p_new[0] << " " << p_new[1] << " " << p_new[2];
    }
  if(!this->_writeBinary)
    {
    output << p[2] << " " << p[1] << " " << p[0];
    }
  else
    {
    this->WriteX<double,float>(output,p[2]);
    this->WriteX<double,float>(output,p[1]);
    this->WriteX<double,float>(output,p[0]);
    }
  ++counter;
}

void VtkWriter::State2Tensor(const State & state, mat_t & D, const int tensorNumber) const
{
  static const size_t local_phi_index = 0;
  static const size_t local_theta_index = 1;
  static const size_t local_psi_index =2;

  vec_t eigenVec1, eigenVec2, eigenVec3;

  if( _full )
    {
    const size_t local_start_index = 6 * (tensorNumber - 1);
    const mat_t & R =
      rotation(
        state[local_start_index + local_phi_index],
        state[local_start_index + local_theta_index],
        state[local_start_index + local_psi_index]);

    // Extract eigenvectors
    eigenVec1 = make_vec(R._[0], R._[3], R._[6]);
    eigenVec2 = make_vec(R._[1], R._[4], R._[7]);
    eigenVec3 = make_vec(R._[2], R._[5], R._[8]);

    }
  else
    {
    // Extract eigenvectors
    eigenVec1 =
      make_vec( state[5 * (tensorNumber - 1) + _p_m1],  state[5 * (tensorNumber - 1) + _p_m2],
                state[5 * (tensorNumber - 1) + _p_m3]);
    eigenVec2 =
      make_vec( state[5 * (tensorNumber - 1) + _p_m1], -state[5 * (tensorNumber - 1) + _p_m2],
                state[5 * (tensorNumber - 1) + _p_m3]);
    eigenVec3 =
      make_vec(-state[5 * (tensorNumber - 1) + _p_m1],  state[5 * (tensorNumber - 1) + _p_m2],
               state[5 * (tensorNumber - 1) + _p_m3]);
    }

  // Perform ijk->RAS transform on eigen vectors
  eigenVec1 = _sizeFreeI2R * eigenVec1;
  eigenVec2 = _sizeFreeI2R * eigenVec2;
  eigenVec3 = _sizeFreeI2R * eigenVec3;

  // Renormalize eigenvectors
  const double vecnorm1 = norm(eigenVec1);
  eigenVec1 = make_vec(eigenVec1._[0] / vecnorm1, eigenVec1._[1] / vecnorm1, eigenVec1._[2] / vecnorm1);
  const double vecnorm2 = norm(eigenVec2);
  eigenVec2 = make_vec(eigenVec2._[0] / vecnorm2, eigenVec2._[1] / vecnorm2, eigenVec2._[2] / vecnorm2);
  const double vecnorm3 = norm(eigenVec3);
  eigenVec3 = make_vec(eigenVec3._[0] / vecnorm3, eigenVec3._[1] / vecnorm3, eigenVec3._[2] / vecnorm3);

  // Compute the diffusion matrix in RAS coordinate system
  // The transformed matrix is still positive-definite
  const mat_t & Q = make_mat(
    eigenVec1._[0], eigenVec2._[0], eigenVec3._[0],
    eigenVec1._[1], eigenVec2._[1], eigenVec3._[1],
    eigenVec1._[2], eigenVec2._[2], eigenVec3._[2]
    );
  D = Q
    * diag(state[5 * (tensorNumber - 1) + _p_l1], state[5 * (tensorNumber - 1) + _p_l2],
           state[5 * (tensorNumber - 1) + _p_l2]) * t(Q) * 1e-6;
}
