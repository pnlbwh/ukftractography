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

#ifndef __qSlicerInteractiveUKFModuleWidget_h
#define __qSlicerInteractiveUKFModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerInteractiveUKFModuleExport.h"

class qSlicerInteractiveUKFModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_INTERACTIVEUKF_EXPORT qSlicerInteractiveUKFModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerInteractiveUKFModuleWidget(QWidget *parent=0);
  virtual ~qSlicerInteractiveUKFModuleWidget();

public slots:


protected:
  QScopedPointer<qSlicerInteractiveUKFModuleWidgetPrivate> d_ptr;

  virtual void setup();

private:
  Q_DECLARE_PRIVATE(qSlicerInteractiveUKFModuleWidget);
  Q_DISABLE_COPY(qSlicerInteractiveUKFModuleWidget);
};

#endif
