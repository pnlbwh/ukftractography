# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName BOOST) #The find_package known name
set(proj        Boost) #This local name
set(${extProjName}_REQUIRED_VERSION "")  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
endif()

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT:STRING=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  ### --- Project specific additions here
  set(Boost_Install_Dir ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install)
  
  if(CMAKE_COMPILER_IS_CLANGXX)
    set(CLANG_ARG -DCMAKE_COMPILER_IS_CLANGXX:BOOL=ON)
  endif()

  ### --- End Project specific additions
  if(CMAKE_COMPILER_IS_CLANGXX)
    set(CLANG_ARG -DCMAKE_COMPILER_IS_CLANGXX:BOOL=ON)
  endif()
  set(BOOST_SOURCE_DIR ${SOURCE_DOWNLOAD_CACHE}/${proj})

  set( Boost_url "http://sourceforge.net/projects/boost/files/boost/1.70.0/boost_1_70_0.tar.gz")
  set( Boost_md5 fea771fe8176828fabf9c09242ee8c26 )
  set( Boost_Bootstrap_Command ./bootstrap.sh )
  set( Boost_b2_Command ./b2 )
  
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(Boost_address_model 64)
else()
  set(Boost_address_model 32)
endif()

 ExternalProject_Add(${proj}
	${${proj}_EP_ARGS}
	BUILD_IN_SOURCE 1
	URL ${Boost_url}
	URL_MD5 ${Boost_md5}
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${Boost_Bootstrap_Command} --prefix=${Boost_Install_Dir}/lib
	BUILD_COMMAND ${Boost_b2_Command} install -j8 --prefix=${Boost_Install_Dir} --with-thread --with-filesystem --with-system --with-date_time --with-program_options  --with-atomic  address-model=${Boost_address_model} link=static,shared 
	INSTALL_COMMAND ""
	)

  set(BOOST_ROOT        ${Boost_Install_Dir})
  set(BOOST_INCLUDE_DIR ${Boost_Install_Dir}/include)
  set(Boost_LIBRARY_DIR ${Boost_Install_Dir}/lib)
  message(STATUS "BOOST_ROOT is " ${BOOST_ROOT})
  message(STATUS "Boost_LIBRARY_DIR is " ${Boost_LIBRARY_DIR})

else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${proj} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  ExternalProject_Add_Empty(${proj} "${${proj}_DEPENDENCIES}")
endif()

mark_as_superbuild(
  VARS
    ${extProjName}_DIR:PATH
  LABELS
    "FIND_PACKAGE"
)

