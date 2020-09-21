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

// .NAME vtkSlicerInteractiveUKFLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerInteractiveUKFLogic_h
#define __vtkSlicerInteractiveUKFLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// STD includes
#include <cstdlib>

#include "vtkSlicerInteractiveUKFModuleLogicExport.h"

class vtkMRMLModelNode;
class vtkMRMLScalarVolumeNode;
class vtkMRMLMarkupsFiducialNode;
class vtkMRMLDiffusionWeightedVolumeNode;
class vtkMRMLCommandLineModuleNode;
class vtkTrivialProducer;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_INTERACTIVEUKF_LOGIC_EXPORT
  vtkSlicerInteractiveUKFLogic : public vtkSlicerModuleLogic
{
public:

  static vtkSlicerInteractiveUKFLogic *New();
  vtkTypeMacro(vtkSlicerInteractiveUKFLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  bool InitTractography(vtkMRMLCommandLineModuleNode*);

  void SetDataNodes(vtkMRMLDiffusionWeightedVolumeNode*,
                    vtkMRMLScalarVolumeNode*,
                    vtkMRMLScalarVolumeNode*,
                    vtkMRMLModelNode*);
  void RunFromSeedPoints(vtkMRMLDiffusionWeightedVolumeNode*,
                         vtkMRMLModelNode*,
                         vtkMRMLMarkupsFiducialNode*);
                         // int pointId);

  void set_seedsPerVoxel(double val);
  void set_stoppingFA(double val);
  void set_seedingThreshold(double val);
  void set_stoppingThreshold(double val);
  void set_numTensor(size_t val);
  void set_stepLength(double val);
  void set_recordLength(double val);
  void set_noddi(bool val);
  void set_freeWater(bool val);

protected:
  vtkSlicerInteractiveUKFLogic();
  virtual ~vtkSlicerInteractiveUKFLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene) override;
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes() override;
  virtual void UpdateFromMRMLScene() override;
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node) override;
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node) override;
private:
  vtkTrivialProducer* producer;

  vtkSlicerInteractiveUKFLogic(const vtkSlicerInteractiveUKFLogic&); // Not implemented
  void operator=(const vtkSlicerInteractiveUKFLogic&); // Not implemented
};

#endif
