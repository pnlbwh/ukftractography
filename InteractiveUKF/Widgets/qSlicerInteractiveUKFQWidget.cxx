/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// Q Widgets includes
#include "qSlicerInteractiveUKFQWidget.h"
#include "ui_qSlicerInteractiveUKFQWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_InteractiveUKF
class qSlicerInteractiveUKFQWidgetPrivate
  : public Ui_qSlicerInteractiveUKFQWidget
{
  Q_DECLARE_PUBLIC(qSlicerInteractiveUKFQWidget);
protected:
  qSlicerInteractiveUKFQWidget* const q_ptr;

public:
  qSlicerInteractiveUKFQWidgetPrivate(
    qSlicerInteractiveUKFQWidget& object);
  virtual void setupUi(qSlicerInteractiveUKFQWidget*);
};

// --------------------------------------------------------------------------
qSlicerInteractiveUKFQWidgetPrivate
::qSlicerInteractiveUKFQWidgetPrivate(
  qSlicerInteractiveUKFQWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerInteractiveUKFQWidgetPrivate
::setupUi(qSlicerInteractiveUKFQWidget* widget)
{
  this->Ui_qSlicerInteractiveUKFQWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerInteractiveUKFQWidget methods

//-----------------------------------------------------------------------------
qSlicerInteractiveUKFQWidget
::qSlicerInteractiveUKFQWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerInteractiveUKFQWidgetPrivate(*this) )
{
  Q_D(qSlicerInteractiveUKFQWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerInteractiveUKFQWidget
::~qSlicerInteractiveUKFQWidget()
{
}
