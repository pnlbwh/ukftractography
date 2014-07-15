get_filename_component(INPUT ${CMAKE_CURRENT_SOURCE_DIR}/../../ukf/Data/Input REALPATH)
get_filename_component(BASELINE ${CMAKE_CURRENT_SOURCE_DIR}/../../ukf/Data/Baseline REALPATH)
set(TESTING_RESULTS_DIRECTORY "${CMAKE_BINARY_DIR}/Testing/Temporary")

file(MAKE_DIRECTORY ${TESTING_RESULTS_DIRECTORY}) ## Need to generate the directory for writing the temporary test outputs to

set(CLP fibertractdispersion)

add_test(NAME fibertractdispersion
  COMMAND $<TARGET_FILE:${CLP}>
  --writeAscii
  --writeUncompressed
  --inputFiberBundle     ${BASELINE}/1T_fiber.vtk
  --outputFiberBundle ${TESTING_RESULTS_DIRECTORY}/${CLP}_1T_fiber.vtk
  )
