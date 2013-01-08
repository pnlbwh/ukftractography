
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
if (NOT ITK_FOUND)
  find_package(ITK REQUIRED)
endif (NOT ITK_FOUND)
include(${ITK_USE_FILE})

if (NOT GenerateCLP_FOUND)
  find_package(GenerateCLP REQUIRED)
endif (NOT GenerateCLP_FOUND)
include(${GenerateCLP_USE_FILE})

# Generate InterfaceCLP.h
SET(vtkFilter_SOURCE vtkFilter.cxx vtkReader.cc vtkWriter.cc FiberFilter.cc Region.cc ExpressionParser.cc ExpressionEvaluator.cc)
GENERATECLP(vtkFilter_SOURCE vtkFilter.xml)

# Add Executables and Link
ADD_EXECUTABLE(vtkFilter ${vtkFilter_SOURCE} )
TARGET_LINK_LIBRARIES(vtkFilter ${ITK_LIBRARIES})
