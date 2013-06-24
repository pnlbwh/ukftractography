cmake_minimum_required(VERSION 2.8.7)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

#-----------------------------------------------------------------------------
set(MODULE_NAME vtkFilter)

#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
#if(NOT Slicer_SOURCE_DIR AND NOT UKFTractography_SuperBuild)
  #find_package(Slicer REQUIRED)
  #include(${Slicer_USE_FILE})
#endif()

#-----------------------------------------------------------------------------
if (UKFTractography_SuperBuild)
  find_package(SlicerExecutionModel REQUIRED)
  include(${SlicerExecutionModel_USE_FILE})

  find_package(ITK REQUIRED)
  include(${ITK_USE_FILE})

  find_package(GenerateCLP REQUIRED)
  include(${GenerateCLP_USE_FILE})
endif()

set(MODULE_SRCS
    vtkFilter.cxx
    vtkWriter.cc
    FiberFilter.cc
    Region.cc
    ExpressionParser.cc
    ExpressionEvaluator.cc
    )

set(MODULE_TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    UKFCommonLib
  )

if(UKFTractography_SuperBuild)
    SEMMacroBuildCLI(
        NAME ${MODULE_NAME}
        TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
        INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
        ADDITIONAL_SRCS ${MODULE_SRCS}
        EXECUTABLE_ONLY
        )
else()
    SEMMacroBuildCLI(
        NAME ${MODULE_NAME}
        TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
        INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
        ADDITIONAL_SRCS ${MODULE_SRCS}
        )
endif()

# Add Executables and Link
#ADD_EXECUTABLE(vtkFilter ${vtkFilter_SOURCE} )
#TARGET_LINK_LIBRARIES(vtkFilter ${ITK_LIBRARIES})
