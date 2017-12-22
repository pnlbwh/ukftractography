/**
 * \file vtk_writer.cc
 * \brief implementation of vtk_writer.h
 * \todo The case differentiation in the beginning is very hackish..
*/

#include "vtkVersion.h"
#include "ukf_types.h"
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
#include "vtkStringArray.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "vtkDataObject.h"

#include "git_version.h"

VtkWriter::VtkWriter(const ISignalData *signal_data, Tractography::model_type filter_model_type, bool write_tensors) :
  _signal_data(signal_data),
  _transform_position(true),
  _filter_model_type(filter_model_type),
  _scale_glyphs(0.01),
  _write_tensors(write_tensors),
  _eigenScaleFactor(1),
  _writeBinary(true),
  _writeCompressed(true)
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
    _p_l1 = 3,
    _p_l2 = 4;
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
  ukfMatrixType i2r = _signal_data->i2r();
  vec3_t              voxel = _signal_data->voxel();
  // Factor out the effect of voxel size
  _sizeFreeI2R <<
    i2r(0, 0) / voxel[2], i2r(0, 1) / voxel[1], i2r(0, 2) / voxel[0],
    i2r(1, 0) / voxel[2], i2r(1, 1) / voxel[1], i2r(1, 2) / voxel[0],
    i2r(2, 0) / voxel[2], i2r(2, 1) / voxel[1], i2r(2, 2) / voxel[0];

}



void VtkWriter
::PopulateFibersAndTensors(vtkPolyData* polyData,
                           const std::vector<UKFFiber>& fibers)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  size_t num_fibers = fibers.size();
  size_t num_points = 0;
  for( size_t i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  for( size_t i = 0; i < num_fibers; ++i )
    {
    size_t fiber_size = fibers[i].position.size();
    for( size_t j = 0; j < fiber_size; ++j )
      {
      vec3_t current = PointConvert(fibers[i].position[j]);
      points->InsertNextPoint(current[0],current[1],current[2]);
      }
    }
  polyData->SetPoints(points);

  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  vtkIdType counter = 0;
  for(size_t i = 0; i < num_fibers; ++i)
    {
    size_t fiber_size = fibers[i].position.size();
    vtkIdType *ids = new vtkIdType[fiber_size];
    for(size_t j = 0; j < fiber_size; ++j)
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
    counter = 0;
    vtkPointData *pointData = polyData->GetPointData();

    mat33_t D;
    for( int local_tensorNumber = 1; local_tensorNumber <= _num_tensors; ++local_tensorNumber )
      {
      vtkSmartPointer<vtkFloatArray> curTensor = vtkSmartPointer<vtkFloatArray>::New();
      curTensor->SetNumberOfComponents(9);
      curTensor->Allocate(num_points * 9);
        {
        std::stringstream ss;
        ss << "tensor" << local_tensorNumber;
        curTensor->SetName(ss.str().c_str());
        }
      for( size_t i = 0; i < num_fibers; i++ )
        {
        const size_t fiber_size = fibers[i].position.size();
        for( size_t j = 0; j < fiber_size; ++j )
          {
          const State & state = fibers[i].state[j];
          State2Tensor(state, D, local_tensorNumber);
          ukfPrecisionType tmp[9];
          for(unsigned ii = 0, v = 0; ii < 3; ++ii)
            {
            for(unsigned jj = 0; jj < 3; ++jj, ++v)
              {
              tmp[v] = D(ii,jj);
              }
            }
          curTensor->InsertNextTuple(tmp);
          }
        }
      const size_t idx = pointData->AddArray(curTensor);
      pointData->SetActiveAttribute(idx,vtkDataSetAttributes::TENSORS);
      }
    }
}

void
VtkWriter
::WritePolyData(vtkSmartPointer<vtkPolyData> pd, const char *filename) const
{
  // Add version information to polydata as vtkStringArray
  std::stringstream header_string;
  header_string << "UKF_GIT_HASH:" << UKF_GIT_HASH;
  vtkSmartPointer<vtkStringArray> version_info = vtkSmartPointer<vtkStringArray>::New();
  version_info->SetNumberOfValues(1);
  version_info->SetValue(0, header_string.str());
  version_info->SetName("UKF_VERSION_INFO");
  pd->GetFieldData()->AddArray(version_info);

  // Output filename extension
  const std::string ext(itksys::SystemTools::GetFilenameExtension(filename));

  if(ext == ".vtp")
    {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    if(this->_writeBinary)
      {
      writer->SetDataModeToBinary();
      }
    else
      {
      writer->SetDataModeToAscii();
      }
    if(this->_writeCompressed)
      {
      writer->SetCompressorTypeToZLib();
      }
    else
      {
      writer->SetCompressorTypeToNone();
      }
#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput(pd);
#else
    writer->SetInputData(pd);
#endif
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
#if (VTK_MAJOR_VERSION < 6)
    writer->SetInput(pd);
#else
    writer->SetInputData(pd);
#endif
    writer->SetFileName(filename);

    writer->Write();
    }
}

int
VtkWriter
::Write(const std::string& file_name,
        const std::string & tractsWithSecondTensor,
        const std::vector<UKFFiber>& fibers,
        bool write_state,
        bool store_glyphs,
        bool if_noddi)
{
  if( fibers.size() == 0 )
    {
    std::cout << "No fiber exists." << std::endl;
    return EXIT_FAILURE;
    }

  //
  // Handle glyps
  if(store_glyphs)
    {
    std::stringstream ss;
    ss << file_name.substr(0, file_name.find_last_of(".") ) << "_glyphs"
       << ".vtk";
    if( WriteGlyphs(ss.str(), fibers)  == EXIT_FAILURE )
      {
      return EXIT_FAILURE;
      }
    }
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  // handle fibers and tensors
  this->PopulateFibersAndTensors(polyData,fibers);

  // no idea if this below is ever called.
  if( !tractsWithSecondTensor.empty() )
    {
    vtkSmartPointer<vtkPolyData> polyData2 = vtkSmartPointer<vtkPolyData>::New();
    this->PopulateFibersAndTensors(polyData2,fibers);
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
  norms->SetName("EstimatedUncertainty");
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
    if(if_noddi)
      fa->SetName("Vic1");
    else
      fa->SetName("FA1");
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
    if(if_noddi)
      fa2->SetName("Vic2");
    else
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
    if(if_noddi)
      trace->SetName("OrientationDispersionIndex1");
    else
      trace->SetName("trace1");
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
    if(if_noddi)
      trace2->SetName("OrientationDispersionIndex2");
    else
      trace2->SetName("trace2");
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
    if(if_noddi)
      free_water->SetName("Viso");
    else
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
    ukfPrecisionType nmse_sum(0);
    unsigned counter(0);
    vtkSmartPointer<vtkFloatArray> normMSE = vtkSmartPointer<vtkFloatArray>::New();
    normMSE->SetNumberOfComponents(1);
    normMSE->Allocate(num_points);
    normMSE->SetName("NormalizedSignalEstimationError");
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
          tmpArray[k] = state[k];
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
  return EXIT_SUCCESS;
}

int VtkWriter::WriteGlyphs(const std::string& file_name,
                            const std::vector<UKFFiber>& fibers)
{
  if( fibers.size() == 0 )
    {
    std::cout << "No fibers existing." << std::endl;
    return EXIT_FAILURE;
    }
  // polyData object to fill in
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int num_fibers = fibers.size();
  int num_points = 0;
  for( int i = 0; i < num_fibers; ++i )
    {
    num_points += fibers[i].position.size();
    }

  int num_tensors = fibers[0].state[0].size() / 5;

  const ukfPrecisionType scale = ukfHalf * _scale_glyphs;


  for( int i = 0; i < num_fibers; ++i )
    {
    int fiber_size = fibers[i].position.size();
    for( int j = 0; j < fiber_size; ++j )
      {
      vec3_t        point = fibers[i].position[j];
      const State& state = fibers[i].state[j];

      // Get the directions.
      vec3_t m1 = vec3_t::Zero();
      vec3_t m2 = vec3_t::Zero();
      vec3_t m3 = vec3_t::Zero();
      if( state.size() == 5 )
        {
        m1 << state[2], state[1], state[0];
        m1 = state[3] / 100.0 * m1;
        }
      else if( state.size() == 6 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        const ukfPrecisionType tmp = m1[0];
        m1[0] = m1[2];
        m1[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        }
      else if( state.size() == 10 )
        {
        m1 << state[2], state[1], state[0];
        m2 << state[7], state[6], state[5];
        m1 = state[3] / 100.0 * m1;
        m2 = state[8] / 100.0 * m2;
        }
      else if( state.size() == 12 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        m2 = rotation_main_dir(state[6], state[7], state[8]);
        ukfPrecisionType tmp = m1[0];
        m1[0] = m1[2];
        m1[2] = tmp;
        tmp = m2[0];
        m2[0] = m2[2];
        m2[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        m2 = state[9] / 100.0 * m2;
        }
      else if( state.size() == 15 )
        {
        m1 << state[2], state[1], state[0];
        m2 << state[7], state[6], state[5];
        m3 <<  state[12], state[11], state[10];
        m1 = state[3] / 100.0 * m1;
        m2 = state[8] / 100.0 * m2;
        m3 = state[13] / 100.0 * m3;
        }
      else if( state.size() == 18 )
        {
        m1 = rotation_main_dir(state[0], state[1], state[2]);
        m2 = rotation_main_dir(state[6], state[7], state[8]);
        m3 = rotation_main_dir(state[12], state[13], state[14]);
        ukfPrecisionType tmp = m1[0];
        m1[0] = m1[2];
        m1[2] = tmp;
        tmp = m2[0];
        m2[0] = m2[2];
        m2[2] = tmp;
        tmp = m3[0];
        m3[0] = m3[2];
        m3[2] = tmp;
        m1 = state[3] / 100.0 * m1;
        m2 = state[9] / 100.0 * m2;
        m3 = state[15] / 100.0 * m3;
        }

      // Calculate the points. The glyphs are represented as two-point lines.
      vec3_t pos1, pos2;
      if( num_tensors == 1 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);
        }
      else if( num_tensors == 2 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);
        }
      else if( num_tensors == 3 )
        {
        pos1 = point - scale * m1;
        pos2 = point + scale * m1;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);

        pos1 = point - scale * m2;
        pos2 = point + scale * m2;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);

        pos1 = point - scale * m3;
        pos2 = point + scale * m3;
        points->InsertNextPoint(pos1[0],pos1[1],pos1[2]);
        points->InsertNextPoint(pos2[0],pos2[1],pos2[2]);
        }
      }
    }
  // do the lines
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

  for( int i = 0; i < num_points * num_tensors; ++i )
    {
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    vtkIdType ids[2];
    ids[0] = i * 2;
    ids[1] = ids[0] + 1;
    line->Initialize(2,ids,points);
    lines->InsertNextCell(line);
    }
  polyData->SetLines(lines);
  WritePolyData(polyData,file_name.c_str());
  return EXIT_SUCCESS;
}

vec3_t
VtkWriter::
PointConvert(const vec3_t& point)
{
  vec3_t rval;
  ukfVectorType p(4);
  p[0] = point[2];    // NOTICE the change of order here. Flips back to the original axis order
  p[1] = point[1];
  p[2] = point[0];

  if( _transform_position )
    {
    p[3] = ukfOne;
    ukfVectorType p_new(4);
    p_new = _signal_data->i2r() * p;    // ijk->RAS transform
    rval[0] = p_new[0];
    rval[1] = p_new[1];
    rval[2] = p_new[2];
    }
  else
    {
    rval[0] = p[2];
    rval[1] = p[1];
    rval[2] = p[0];
    }
  return rval;
}

void VtkWriter::State2Tensor(const State & state, mat33_t & D, const int tensorNumber) const
{
  vec3_t eigenVec1;
  const int tensorIndex= ( tensorNumber - 1 );
  const size_t after_tensor_offset_index = _tensor_space * (tensorIndex);


  if( _full )
    {
    static const size_t local_phi_index = 0;
    static const size_t local_theta_index = 1;
    static const size_t local_psi_index =2;

    const mat33_t & R =
      rotation(
        state[after_tensor_offset_index + local_phi_index],
        state[after_tensor_offset_index + local_theta_index],
        state[after_tensor_offset_index + local_psi_index]);

    // Extract eigenvectors
    eigenVec1 << R(0,0), R(1,0), R(2,0);
/* HACK THIS SEEMS WRONG!  Why have eigenVec2 or eigenVec3?
    vec3_t eigenVec2;
    vec3_t eigenVec3;
    eigenVec2 << R(0,1), R(1,1), R(2,1);
    eigenVec3 << R(0,2), R(1,2), R(2,2);
*/
    }
  else
  {
      eigenVec1 <<
          state[after_tensor_offset_index + _p_m1],
          state[after_tensor_offset_index + _p_m2],
          state[after_tensor_offset_index + _p_m3];

      // Perform ijk->RAS transform on eigen vectors
      eigenVec1 = _sizeFreeI2R * eigenVec1;

      // Renormalize eigenvectors
      eigenVec1.normalize();
  }

  // Compute the diffusion matrix in RAS coordinate system
  // The transformed matrix is still positive-definite
  const float L1= state[after_tensor_offset_index + _p_l1];
  const float L2= state[after_tensor_offset_index + _p_l2];
#if 1
  D = diffusion_l2eql3(eigenVec1, L1, L2);
#else
  diagmat3_t lambdas;
  lambdas.diagonal()[0] = L1;
  lambdas.diagonal()[1] = L2;
  lambdas.diagonal()[2] = lambdas.diagonal()[1];
  D=diffusion(eigenVec1,lambdas);
#endif
}
