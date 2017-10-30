from __future__ import print_function

import sys, os

#
# EVERYTHING BELOW THIS LINE MUST BE IN TRY BLOCK
#

# wrap code in try block so that we can exit with error if something
# unexpectedly goes wrong. otherwise Slicer may sit open until
# killed by CTest timeout.
try:
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
    sys.exit(1)

except Exception as err:
  print("Error: caught Python exception!: ", err)
  sys.exit(1)

#
# EVERYTHING ABOVE THIS LINE MUST BE IN TRY BLOCK
#

# clean completion, return EXIT_SUCCESS
# note use Python sys.exit because slicer.util.exit is effectively
# asynchronous: https://issues.slicer.org/view.php?id=4470
sys.exit(0)
