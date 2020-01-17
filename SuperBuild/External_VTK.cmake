
set(proj VTK)

set(VTK_VERSION_MAJOR 7)

# Set dependency list
set(${proj}_DEPENDENCIES "zlib")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(VTK_DIR CACHE)
  unset(VTK_SOURCE_DIR CACHE)
  find_package(VTK REQUIRED NO_MODULE)
endif()

# Sanity checks
if(DEFINED VTK_DIR AND NOT EXISTS ${VTK_DIR})
  message(FATAL_ERROR "VTK_DIR variable is defined but corresponds to non-existing directory")
endif()

if(DEFINED VTK_SOURCE_DIR AND NOT EXISTS ${VTK_SOURCE_DIR})
  message(FATAL_ERROR "VTK_SOURCE_DIR variable is defined but corresponds to non-existing directory")
endif()

# For MinGW for case of compilation failure cause of 'too many sections' error
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  if(MINGW)
	set(ep_common_cxx_flags "${ep_common_cxx_flags} -Wa,-mbig-obj")
  elseif(MSVC)
	set(ep_common_cxx_flags "${ep_common_cxx_flags} /bigobj")
  endif()
endif()

if((NOT DEFINED VTK_DIR OR NOT DEFINED VTK_SOURCE_DIR) AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  set(${proj}_GIT_REPOSITORY "${git_protocol}://github.com/Kitware/VTK.git" CACHE STRING "Repository from which to get VTK" FORCE)
  set(${proj}_GIT_TAG "v7.1.1") #"b86da7eef93f75c4a7f524b3644523ae6b651bc4")  # VTK v7.1.1

## Use ../VTK/Utilities/Maintenance/WhatModulesVTK.py ../VTK ./
## to identify necessary modules for VTK

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj}
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build
    GIT_REPOSITORY "${${proj}_GIT_REPOSITORY}"
    GIT_TAG ${${proj}_GIT_TAG}
    CMAKE_ARGS -Wno-dev --no-warn-unused-cli
    CMAKE_CACHE_ARGS
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DCMAKE_INCLUDE_DIRECTORIES_BEFORE:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DVTK_WRAP_TCL:BOOL=OFF
      -DVTK_WRAP_PYTHON:BOOL=OFF
      -DVTK_USE_GUISUPPORT:BOOL=OFF
      -DVTK_USE_QT:BOOL=OFF
      -DVTK_BUILD_ALL_MODULES_FOR_TESTS:BOOL=OFF
      -DVTK_Group_Rendering:BOOL=OFF
      -DVTK_Group_StandAlone:BOOL=OFF
      -DModule_vtkCommonCore:BOOL=ON
      -DModule_vtkCommonDataModel:BOOL=ON
      -DModule_vtkCommonExecutionModel:BOOL=ON
      -DModule_vtkCommonMath:BOOL=ON
      -DModule_vtkCommonMisc:BOOL=ON
      -DModule_vtkCommonSystem:BOOL=ON
      -DModule_vtkCommonTransforms:BOOL=ON
      -DModule_vtkIOCore:BOOL=ON
      -DModule_vtkIOGeometry:BOOL=ON
      -DModule_vtkIOLegacy:BOOL=ON
      -DModule_vtkIOXML:BOOL=ON
      -DModule_vtkIOXMLParser:BOOL=ON
      INSTALL_COMMAND ""
    )
  ### --- End Project specific additions
  set(${proj}_DIR ${CMAKE_CURRENT_BINARY_DIR}/${proj}-build)

  set(VTK_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})

  set(PNG_INCLUDE_DIR ${VTK_SOURCE_DIR}/Utilities/vtkpng)

  set(PNG_LIBRARY_DIR ${VTK_DIR}/bin)
  if(CMAKE_CONFIGURATION_TYPES)
    set(PNG_LIBRARY_DIR ${PNG_LIBRARY_DIR}/${CMAKE_CFG_INTDIR})
  endif()
  if(WIN32)
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/vtkpng.lib)
  elseif(APPLE)
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/libvtkpng.dylib)
  else()
    set(PNG_LIBRARY ${PNG_LIBRARY_DIR}/libvtkpng.so)
  endif()

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDENCIES})
endif()

mark_as_superbuild(VTK_SOURCE_DIR:PATH)

mark_as_superbuild(
  VARS ${proj}_DIR:PATH VTK_VERSION_MAJOR:STRING
  LABELS "FIND_PACKAGE"
  )
