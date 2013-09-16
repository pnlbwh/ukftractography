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
#include "itksys/SystemTools.hxx"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkLine.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkDataObject.h"
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
              this->WriteX<double,float>(output,D._[k] * _eigenScaleFactor);
              }
            }
          }
        }
      }
    }
  output << std::endl;
}

void VtkWriter
::PopulateFibersAndTensors(vtkSmartPointer<vtkPolyData> &polyData,
                           const std::vector<UKFFiber>& fibers)
{
  polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j )
      {
      vec_t current = PointConvert(fibers[i].position[j]);
      points->InsertNextPoint(current._);
      }
    }
  polyData->SetPoints(points);

  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  vtkIdType counter = 0;
  for(int i = 0; i < num_fibers; ++i)
    {
    int fiber_size = fibers[i].position.size();
    vtkIdType *ids = new vtkIdType[fiber_size];
    for(int j = 0; j < fiber_size; ++j)
      {
      ids[j] = counter;
      counter++;
      }
    vtkSmartPointer<vtkLine> curLine = vtkSmartPointer<vtkLine>::New();
    curLine->Initialize(fiber_size,ids,points);
    lines->InsertNextCell(curLine);
    delete [] ids;
    }
  polyData->SetLines(lines);

  /////
  // Dataset attribute part starts
  /////


  if( _write_tensors )
    {
    typedef std::vector<double> State;
    counter = 0;
    vtkPointData *pointData = polyData->GetPointData();

    for( int local_tensorNumber = 1; local_tensorNumber <= _num_tensors; ++local_tensorNumber )
      {
      vtkSmartPointer<vtkFloatArray> curTensor = vtkSmartPointer<vtkFloatArray>::New();
      curTensor->SetNumberOfComponents(9);
      curTensor->Allocate(num_points * 9);
      std::stringstream ss;
      ss << local_tensorNumber;
      curTensor->SetName(ss.str().c_str());
      // output << "TENSORS Tensor" << local_tensorNumber << " float" << std::endl;
      for( int i = 0; i < num_fibers; i++ )
        {
        const int fiber_size = static_cast<int>(fibers[i].position.size() );
        for( int j = 0; j < fiber_size; ++j )
          {
          const State & state = fibers[i].state[j];
          mat_t         D;
          State2Tensor(state, D, local_tensorNumber);
          curTensor->InsertNextTuple(D._);
          }
        }
      int idx = pointData->AddArray(curTensor);
      pointData->SetActiveAttribute(idx,vtkDataSetAttributes::TENSORS);
      }
    }
}

void
VtkWriter
::WritePolyData(const vtkPolyData *pd, const char *filename) const
{
  const std::string ext(itksys::SystemTools::GetFilenameExtension(filename));
  // all the casting is because vtk SetInput isn't const-correct.
  vtkDataObject *dataObject =
    const_cast<vtkDataObject *>(static_cast<const vtkDataObject *>(pd));
  if(ext == ".vtp")
    {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetDataModeToBinary();
    writer->SetCompressorTypeToZLib();
    writer->SetInput(dataObject);
    writer->SetFileName(filename);
    writer->Write();
    }
  else
    {
    vtkSmartPointer<vtkPolyDataWriter> writer =
      vtkSmartPointer<vtkPolyDataWriter>::New();
    if(this->_writeBinary)
      {
      writer->SetFileTypeToBinary();
      }
    writer->SetInput(dataObject);
    writer->SetFileName(filename);
    writer->Write();
    }
}

bool
VtkWriter::
Write(const std::string& file_name,
      const std::string & tractsWithSecondTensor,
      const std::vector<UKFFiber>& fibers,
      bool write_state,
      bool store_glyphs)
{
  if( fibers.size() == 0 )
    {
    std::cout << "No fiber exists." << std::endl;
    return true;
    }

  //
  // Handle glyps
  if(store_glyphs)
    {
    std::stringstream ss;
    ss << file_name.substr(0, file_name.find_last_of(".") ) << "_glyphs"
       << ".vtk";
    if( WriteGlyphs(ss.str(), fibers) )
      {
      return true;
      }
    }
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData;
  // handle fibers and tensors
  this->PopulateFibersAndTensors(polyData,fibers);

  // no idea if this below is ever called.
  if( !tractsWithSecondTensor.empty() )
    {
    vtkSmartPointer<vtkPolyData> polyData2;
    this->PopulateFibersAndTensors(polyData,fibers);
    WritePolyData(polyData2, tractsWithSecondTensor.c_str());
    }

  // norm, fa etc hung as arrays on the point data for the polyData
  // object.
  vtkPointData *pointData = polyData->GetPointData();

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  // write norm
  {
  vtkSmartPointer<vtkFloatArray> norms = vtkSmartPointer<vtkFloatArray>::New();
  norms->SetNumberOfComponents(1);
  norms->Allocate(num_points);
  norms->SetName("norm");
  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j)
      {
      norms->InsertNextValue(fibers[i].norm[j]);
      }
    }
  int idx = pointData->AddArray(norms);
  pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
  }

  // write fa
  if(fibers[0].fa.size() > 0)
    {
    vtkSmartPointer<vtkFloatArray> fa = vtkSmartPointer<vtkFloatArray>::New();
    fa->SetNumberOfComponents(1);
    fa->Allocate(num_points);
    fa->SetName("FA");
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        fa->InsertNextValue(fibers[i].fa[j]);
        }
      }
    int idx = pointData->AddArray(fa);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    }

  // fa2
  if(fibers[0].fa2.size() > 0)
    {
    vtkSmartPointer<vtkFloatArray> fa2 = vtkSmartPointer<vtkFloatArray>::New();
    fa2->SetName("FA2");
    fa2->SetNumberOfComponents(1);
    fa2->Allocate(num_points);
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        fa2->InsertNextValue(fibers[i].fa2[j]);
        }
      }
    int idx = pointData->AddArray(fa2);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    }

  // trace
  if(fibers[0].trace.size() > 0)
    {
    vtkSmartPointer<vtkFloatArray> trace = vtkSmartPointer<vtkFloatArray>::New();
    trace->SetNumberOfComponents(1);
    trace->Allocate(num_points);
    trace->SetName("Trace");
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        trace->InsertNextValue(fibers[i].trace[j]);
        }
      }
    int idx = pointData->AddArray(trace);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    }

  // trace2
  if(fibers[0].trace2.size() > 0)
    {
    vtkSmartPointer<vtkFloatArray> trace2 = vtkSmartPointer<vtkFloatArray>::New();
    trace2->SetNumberOfComponents(1);
    trace2->Allocate(num_points);
    trace2->SetName("Trace2");
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        trace2->InsertNextValue(fibers[i].trace2[j]);
        }
      }
    int idx = pointData->AddArray(trace2);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    }

  if(fibers[0].free_water.size() > 0)
    {
    vtkSmartPointer<vtkFloatArray> free_water = vtkSmartPointer<vtkFloatArray>::New();
    free_water->SetNumberOfComponents(1);
    free_water->Allocate(num_points);
    free_water->SetName("FreeWater");
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        free_water->InsertNextValue(fibers[i].free_water[j]);
        }
      }
    int idx = pointData->AddArray(free_water);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    }

  if(fibers[0].normMSE.size() > 0)
    {
    double nmse_sum(0);
    unsigned counter(0);
    vtkSmartPointer<vtkFloatArray> normMSE = vtkSmartPointer<vtkFloatArray>::New();
    normMSE->SetNumberOfComponents(1);
    normMSE->Allocate(num_points);
    normMSE->SetName("NMSE");
    for( int i = 0; i < num_fibers; ++i )
      {
      int fiber_size = fibers[i].position.size();
      for( int j = 0; j < fiber_size; ++j )
        {
        normMSE->InsertNextValue(fibers[i].normMSE[j]);
        nmse_sum += fibers[i].normMSE[j];
        ++counter;
        }
      }
    int idx = pointData->AddArray(normMSE);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::SCALARS);
    std::cout << "nmse_avg=" << nmse_sum / counter << std::endl;
    }
  else
    {
    std::cout << "nmse_avg=0" << std::endl;
    }

  if(write_state)
    {
    int state_dim = fibers[0].state[0].size();
    vtkSmartPointer<vtkFloatArray> stateArray = vtkSmartPointer<vtkFloatArray>::New();
    stateArray->SetNumberOfComponents(state_dim);
    stateArray->Allocate(num_points * state_dim);
    stateArray->SetName("state");

    float *tmpArray = new float[state_dim];

    for( int i = 0; i < num_fibers; i++ )
      {
      const int fiber_size = static_cast<int>(fibers[i].position.size() );
      for( int j = 0; j < fiber_size; ++j )
        {
        const State & state = fibers[i].state[j];
        for(int k = 0; k < state_dim; ++k)
          {
          tmpArray[k] = state[i];
          }
        stateArray->InsertNextTuple(tmpArray);
        }
      }
    int idx = pointData->AddArray(stateArray);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::VECTORS);
    delete [] tmpArray;
    }

  if(fibers[0].covariance.size() > 0)
    {
    int state_dim = fibers[0].state[0].size();

    int cov_dim = state_dim * (state_dim + 1) / 2;

    vtkSmartPointer<vtkFloatArray> covarianceArray = vtkSmartPointer<vtkFloatArray>::New();
    covarianceArray->SetNumberOfComponents(cov_dim);
    covarianceArray->Allocate(num_points * cov_dim);
    covarianceArray->SetName("covariance");

    float *tmpArray = new float[cov_dim];

    for( int i = 0; i < num_fibers; i++ )
      {
      const int fiber_size = static_cast<int>(fibers[i].position.size() );
      for( int j = 0; j < fiber_size; ++j )
        {
        int covIndex = 0;
        for(int a = 0; a < state_dim; ++a)
          {
          for(int b = a; b < state_dim; ++b)
            {
            tmpArray[covIndex] = fibers[i].covariance[j](a,b);
            ++covIndex;
            }
          }
        covarianceArray->InsertNextTuple(tmpArray);
        }
      }
    int idx = pointData->AddArray(covarianceArray);
    pointData->SetActiveAttribute(idx,vtkDataSetAttributes::VECTORS);
    delete [] tmpArray;
    }

  WritePolyData(polyData,file_name.c_str());
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
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  int num_tensors = fibers[0].state[0].size() / 5;

  const double scale = 0.5 * _scale_glyphs;


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
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);
        }
      else if( num_tensors == 2 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);
        }
      else if( num_tensors == 3 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);

        pos1 = point - scale * m3;
        pos2 = point + scale * m3;
        points->InsertNextPoint(pos1._);
        points->InsertNextPoint(pos2._);
        }
      }
    }
  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  for( int i = 0; i < num_points * num_tensors; ++i )
    {
    // output << "2 " << i * 2 << " " << i * 2 + 1 << std::endl;
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    vtkIdType ids[2];
    ids[0] = i * 2;
    ids[1] = ids[0] + 1;
    line->Initialize(2,ids,points);
    lines->InsertNextCell(line);
    }
  polyData->SetLines(lines);
  WritePolyData(polyData,file_name.c_str());
  return false;
}

vec_t
VtkWriter::
PointConvert(const vec_t& point)
{
  vec_t rval;
  vnl_vector<double> p(4);
  p[0] = point._[2];    // NOTICE the change of order here. Flips back to the original axis order
  p[1] = point._[1];
  p[2] = point._[0];

  if( _transform_position )
    {
    p[3] = 1.0;
    vnl_vector<double> p_new(4);
    p_new = _signal_data->i2r() * p;    // ijk->RAS transform
    rval._[0] = p_new[0];
    rval._[1] = p_new[1];
    rval._[2] = p_new[2];
    }
  else
    {
    rval._[0] = p[2];
    rval._[1] = p[1];
    rval._[2] = p[0];
    }
  return rval;
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
