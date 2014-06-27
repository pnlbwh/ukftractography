
#-----------------------------------------------------------------------------
set(MODULE_NAME vtkFilter)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

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

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  )
