from __future__ import print_function

import sys, os


#
# Get widgets and set test parameters
#

cli_widget = slicer.modules.ukftractography.widgetRepresentation()

dwi_selector = findChild(cli_widget, "dwiFile")
tck_selector = findChild(cli_widget, "tracts")
msk_selector = findChild(cli_widget, "maskFile")
lbl_selector = findChild(cli_widget, "seedsFile")

# make sure selectors exist
if not (tck_selector and dwi_selector and tck_selector and lbl_selector):
  print("Error: Unable to get selector combo box(es) from CLI widget!", file=sys.stderr)
  slicer.util.exit(1)

#
# Load test data
#

slicer.util.loadVolume(os.path.join(data_dir, "two_tensor_fw.nhdr"))
slicer.util.loadVolume(os.path.join(data_dir, "mask.nhdr"))
slicer.util.loadVolume(os.path.join(data_dir, "seed.nhdr"))

dwi_selector.setCurrentNode(getNode("two_tensor_fw"))
msk_selector.setCurrentNode(getNode("mask"))
lbl_selector.setCurrentNode(getNode("seed"))

# Hack the node type here and create output node.
# Avoids dependency for vtkMRMLFiberBundleNode.

tck_selector.nodeTypes = (u'vtkMRMLModelNode',)
tck_selector.addNode()

#
# Run the CLI
#
cli_mod_node = cli_widget.currentCommandLineModuleNode()
# run synchronously
cli_widget.module().cliModuleLogic().ApplyAndWait(cli_mod_node)
result = cli_mod_node.GetStatus()

#
# Return exit status for CTest
#

if result != cli_mod_node.Completed:
  print("Error: CLI completed with errors!", file=sys.stderr)
  slicer.util.exit(1)

slicer.util.exit(0)
