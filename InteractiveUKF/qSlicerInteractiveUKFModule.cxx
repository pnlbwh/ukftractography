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

// Qt includes
#include <QtPlugin>

// InteractiveUKF Logic includes
#include <vtkSlicerInteractiveUKFLogic.h>

// InteractiveUKF includes
#include "qSlicerInteractiveUKFModule.h"
#include "qSlicerInteractiveUKFModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerInteractiveUKFModule, qSlicerInteractiveUKFModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerInteractiveUKFModulePrivate
{
public:
  qSlicerInteractiveUKFModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerInteractiveUKFModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerInteractiveUKFModulePrivate::qSlicerInteractiveUKFModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerInteractiveUKFModule methods

//-----------------------------------------------------------------------------
qSlicerInteractiveUKFModule::qSlicerInteractiveUKFModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerInteractiveUKFModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerInteractiveUKFModule::~qSlicerInteractiveUKFModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerInteractiveUKFModule::helpText() const
{
  return "This is a loadable module that can be bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerInteractiveUKFModule::acknowledgementText() const
{
  return "This work was partially funded by NIH grant NXNNXXNNNNNN-NNXN";
}

//-----------------------------------------------------------------------------
QStringList qSlicerInteractiveUKFModule::contributors() const
{
  QStringList moduleContributors;
  moduleContributors << QString("John Doe (AnyWare Corp.)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerInteractiveUKFModule::icon() const
{
  return QIcon(":/Icons/InteractiveUKF.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerInteractiveUKFModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerInteractiveUKFModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerInteractiveUKFModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation* qSlicerInteractiveUKFModule
::createWidgetRepresentation()
{
  return new qSlicerInteractiveUKFModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerInteractiveUKFModule::createLogic()
{
  return vtkSlicerInteractiveUKFLogic::New();
}
