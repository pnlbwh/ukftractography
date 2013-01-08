
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
message("exec output path: ${EXECUTABLE_OUTPUT_PATH}")

if (NOT ITK_FOUND)
  find_package(ITK REQUIRED)
endif (NOT ITK_FOUND)
include(${ITK_USE_FILE})

if (NOT GenerateCLP_FOUND)
  find_package(GenerateCLP REQUIRED)
endif (NOT GenerateCLP_FOUND)
include(${GenerateCLP_USE_FILE})

# Generate InterfaceCLP.h
SET(vtk2mask_SOURCE vtk2mask.cxx vtkReader.cc Converter.cc) 
GENERATECLP(vtk2mask_SOURCE vtk2mask.xml)

# Add Executables and Link
ADD_EXECUTABLE(vtk2mask ${vtk2mask_SOURCE} )
TARGET_LINK_LIBRARIES(vtk2mask ${ITK_LIBRARIES})

