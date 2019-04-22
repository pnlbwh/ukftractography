/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// InteractiveUKF Logic includes
#include "vtkSlicerInteractiveUKFLogic.h"

// SEM includes
#include <ModuleDescription.h>

// Slicer includes
#include <vtkMRMLScene.h>
#include <vtkTeemNRRDWriter.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include <vtkMRMLDiffusionWeightedVolumeNode.h>
#include <vtkMRMLCommandLineModuleNode.h>
#include <vtkSlicerCLIModuleLogic.h>

// ITK includes
#include <itkVersion.h>
#if ITK_VERSION_MAJOR >= 5
#include <itkMultiThreaderBase.h>
#else
#include <itkMultiThreader.h>
#endif

// VTK includes
#include <vtkNew.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkAlgorithm.h>
#include <vtkAlgorithmOutput.h>
#include <vtkTrivialProducer.h>

// Teem includes
#include "teem/nrrd.h"

// STD includes
#include <cassert>

// UKF includes
#include <tractography.h>
#include <cli.h>
#include "BRAINSThreadControl.h"

static Tractography* g_tracto;

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerInteractiveUKFLogic);

//----------------------------------------------------------------------------
vtkSlicerInteractiveUKFLogic::vtkSlicerInteractiveUKFLogic() :
  producer(vtkTrivialProducer::New())
{
}

//----------------------------------------------------------------------------
vtkSlicerInteractiveUKFLogic::~vtkSlicerInteractiveUKFLogic()
{
  producer->Delete();
}

//----------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

static int CLILoader(int argc, char** argv)
{
  UKFSettings ukf_settings;
  ukf_parse_cli(argc, argv, ukf_settings);

  // NOTE:  When used as share libary one must be careful not to permanently reset number of threads
  //        for entire program (i.e. when used as a slicer modules.
  //        This also addresses the issue when the program is run as part of a batch processing
  //        system so that the number of cores allocated by scheduler is respected rather than
  //        blindly using all the cores that are found.
  //        This implementation is taken from extensive testing of the BRAINSTools
  // (this object will be deleted by RAII and return to original thread count)
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(ukf_settings.num_threads);
#if ITK_VERSION_MAJOR >= 5
  const int actualNumThreadsUsed = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();
#else
  const int actualNumThreadsUsed = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
#endif
  ukf_settings.num_threads = actualNumThreadsUsed;
  {
    std::cout << "Found " << actualNumThreadsUsed << " cores on your system." << std::endl;
    std::cout << "Running tractography with " << actualNumThreadsUsed << " thread(s)." << std::endl;
  }

  g_tracto = new Tractography(ukf_settings);

  return EXIT_SUCCESS;
}

//---------------------------------------------------------------------------
bool vtkSlicerInteractiveUKFLogic::InitTractography(vtkMRMLCommandLineModuleNode* n)
{
  assert(this->GetMRMLScene());
  assert(this->GetMRMLApplicationLogic());

  std::stringstream entrypoint;
  entrypoint << "slicer:";
  entrypoint << reinterpret_cast<void*>(&CLILoader);

  vtkNew<vtkMRMLCommandLineModuleNode> cli;
  cli->Copy(n);
  cli->GetModuleDescription().SetTarget(entrypoint.str());
  cli->GetModuleDescription().SetType("SharedObjectModule");

  /* Make a copy of the node to run directly */
  vtkNew<vtkSlicerCLIModuleLogic> cli_logic;
  cli_logic->SetMRMLScene(this->GetMRMLScene());
  cli_logic->SetMRMLApplicationLogic(this->GetMRMLApplicationLogic());
  cli_logic->ApplyAndWait(cli.GetPointer(), false);

  return EXIT_SUCCESS;
}
//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void setWriterProps(vtkMRMLVolumeNode* vol, vtkTeemNRRDWriter* writer)
{
  vtkNew<vtkMatrix4x4> ijkToRas;
  vol->GetIJKToRASMatrix(ijkToRas.GetPointer());
  writer->SetIJKToRASMatrix(ijkToRas.GetPointer());
}


//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::SetDataNodes(
    vtkMRMLDiffusionWeightedVolumeNode* dwiNode,
    vtkMRMLScalarVolumeNode* maskNode,
    vtkMRMLScalarVolumeNode* seedNode,
    vtkMRMLModelNode* fbNode)
{
  assert(dwiNode->IsA("vtkMRMLDiffusionWeightedVolumeNode"));
  assert(maskNode->IsA("vtkMRMLScalarVolumeNode"));
  assert(fbNode->IsA("vtkMRMLFiberBundleNode"));

  vtkNew<vtkTeemNRRDWriter> writer;

  writer->SetInputConnection(dwiNode->GetImageDataConnection());
  setWriterProps((vtkMRMLVolumeNode*)dwiNode, writer.GetPointer());
  vtkNew<vtkMatrix4x4> mf;
  dwiNode->GetMeasurementFrameMatrix(mf.GetPointer());
  writer->SetMeasurementFrameMatrix(mf.GetPointer());

  vtkDoubleArray* grads = NULL;
  vtkDoubleArray* bValues = NULL;
  grads = dwiNode->GetDiffusionGradients();
  bValues = dwiNode->GetBValues();

  if (grads)
    writer->SetDiffusionGradients(grads);
  if (bValues)
    writer->SetBValues(bValues);

  Nrrd* nrrd = (Nrrd*)writer->MakeNRRD();

  vtkNew<vtkTeemNRRDWriter> maskWriter;
  maskWriter->SetInputConnection(maskNode->GetImageDataConnection());
  setWriterProps((vtkMRMLVolumeNode*)maskNode, maskWriter.GetPointer());

  Nrrd* mask = (Nrrd*)maskWriter->MakeNRRD();

  vtkNew<vtkTeemNRRDWriter> seedWriter;
  Nrrd* seed = NULL;
  if (seedNode)
    {
    seedWriter->SetInputConnection(seedNode->GetImageDataConnection());
    setWriterProps((vtkMRMLVolumeNode*)seedNode, seedWriter.GetPointer());
    seed = (Nrrd*)seedWriter->MakeNRRD();
    }

  Tractography* tract = dynamic_cast<Tractography*>(g_tracto);
  if (!tract)
    {
    std::cerr << "Uninitialized UKFTractography!" << std::endl;
    return;
    }

  tract->SetData(nrrd, mask, seed, false /*normalizedDWIData*/);
  tract->UpdateFilterModelType();

  vtkPolyData* pd = vtkPolyData::New();
  //SafeDownCast(producer->GetOutputDataObject(0));
  this->producer->SetOutput(pd);
  pd->Delete();

  fbNode->SetPolyDataConnection(this->producer->GetOutputPort());
}

//---------------------------------------------------------------------------
void vtkSlicerInteractiveUKFLogic::RunFromSeedPoints
      (vtkMRMLDiffusionWeightedVolumeNode* dwiNode,
       vtkMRMLModelNode* fbNode,
       vtkMRMLMarkupsFiducialNode* markupsNode)
{
#ifdef NDEBUG
  (void)fbNode;
#else
  assert(fbNode->IsA("vtkMRMLFiberBundleNode"));
#endif
  assert(markupsNode->IsA("vtkMRMLMarkupsNode"));
  assert(dwiNode->IsA("vtkMRMLDiffusionWeightedVolumeNode"));

  Tractography* tract = dynamic_cast<Tractography*>(g_tracto);
  if (!tract)
    {
    std::cerr << "Uninitialized UKFTractography!" << std::endl;
    return;
    }

  vtkPolyData* pd = vtkPolyData::SafeDownCast(producer->GetOutputDataObject(0));
  if (!pd)
    {
    std::cerr << "No polydata!" << std::endl;
    assert(pd);
    return;
    }
  // avoid crashes and sync issues during update
  pd->Reset();

  vtkNew<vtkTransform> RASxfmIJK;
  vtkNew<vtkMatrix4x4> RAStoIJK;
  dwiNode->GetRASToIJKMatrix(RAStoIJK.GetPointer());
  RASxfmIJK->SetMatrix(RAStoIJK.GetPointer());
  stdVec_t seeds;

  for (int i = 0; i < markupsNode->GetNumberOfFiducials(); i++)
    {
    vec3_t pos_in, pos_out;
    markupsNode->GetNthFiducialPosition(i, pos_in.data());
    RASxfmIJK->TransformPoint(pos_in.data(), pos_out.data());
    pos_out = vec3_t(pos_out[2], pos_out[1], pos_out[0]); // axis order for nrrd
    seeds.push_back(pos_out);
    }

  tract->SetSeeds(seeds);

  tract->SetOutputPolyData(pd);
  tract->Run();

  // TODO fix
  // work around https://issues.slicer.org/view.php?id=3786
  this->producer->Update();
  this->producer->Modified();
}

void vtkSlicerInteractiveUKFLogic::set_seedsPerVoxel(double val) {
  if (!g_tracto) return;
  g_tracto->_seeds_per_voxel = val;
}

void vtkSlicerInteractiveUKFLogic::set_stoppingFA(double val){
  if (!g_tracto) return;
  g_tracto->_fa_min = val;
}
void vtkSlicerInteractiveUKFLogic::set_seedingThreshold(double val) {
  if (!g_tracto) return;
  g_tracto->_seeding_threshold = val;
}
void vtkSlicerInteractiveUKFLogic::set_stoppingThreshold(double val) {
  if (!g_tracto) return;
  g_tracto->_mean_signal_min = val;
}
void vtkSlicerInteractiveUKFLogic::set_stepLength(double val) {
  if (!g_tracto) return;
  g_tracto->_stepLength = val;
}
void vtkSlicerInteractiveUKFLogic::set_recordLength(double val) {
  if (!g_tracto) return;
  g_tracto->_steps_per_record = val/g_tracto->_stepLength;
}

// these calls invalidate the model and must call UpdateFilterModelType()
void vtkSlicerInteractiveUKFLogic::set_numTensor(size_t val) {
  if (!g_tracto) return;
  g_tracto->_num_tensors = val;
  g_tracto->UpdateFilterModelType();
}
void vtkSlicerInteractiveUKFLogic::set_noddi(bool val) {
  if (!g_tracto) return;
  g_tracto->_noddi = val;
  g_tracto->UpdateFilterModelType();
}

void vtkSlicerInteractiveUKFLogic::set_freeWater(bool val) {
  if (!g_tracto) return;
  g_tracto->_free_water = val;
  g_tracto->UpdateFilterModelType();
}
