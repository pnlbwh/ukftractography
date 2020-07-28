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

  if(UNIX)
    set(Boost_url "http://sourceforge.net/projects/boost/files/boost/1.70.0/boost_1_70_0.tar.gz")
    set(Boost_md5 fea771fe8176828fabf9c09242ee8c26)
    set(Boost_Bootstrap_Command ./bootstrap.sh)
    set(Boost_b2_Command ./b2)
  else()
    if(WIN32)
      set(Boost_url "http://sourceforge.net/projects/boost/files/boost/1.70.0/boost_1_70_0.zip")
      set(Boost_md5 a110ebd91a3d2c34c72ace09c92ae50b)
      set(Boost_Bootstrap_Command ./bootstrap.bat)
      set(Boost_b2_Command b2)
    endif()
  endif()

  if(MSVC)
    if(MSVC_VERSION GREATER_EQUAL 1400 AND MSVC_VERSION LESS 1500)
	  list(APPEND Boost_b2_Command toolset=msvc-8.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1500 AND MSVC_VERSION LESS 1600)
	  list(APPEND Boost_b2_Command toolset=msvc-9.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1600 AND MSVC_VERSION LESS 1700)
	  list(APPEND Boost_b2_Command toolset=msvc-10.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1700 AND MSVC_VERSION LESS 1800)
	  list(APPEND Boost_b2_Command toolset=msvc-11.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1800 AND MSVC_VERSION LESS 1900)
	  list(APPEND Boost_b2_Command toolset=msvc-12.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1900 AND MSVC_VERSION LESS 1910)
	  list(APPEND Boost_b2_Command toolset=msvc-14.0)
    elseif(MSVC_VERSION GREATER_EQUAL 1910 AND MSVC_VERSION LESS 1920)
	  list(APPEND Boost_b2_Command toolset=msvc-14.1)
    elseif(MSVC_VERSION GREATER_EQUAL 1920 AND MSVC_VERSION LESS 1927)
	  list(APPEND Boost_b2_Command toolset=msvc-14.2)
    else()	
	  message(FATAL_ERROR "Unknown MSVC compiler version [${MSVC_VERSION}]")
	endif()
  endif()

  if(XCODE_VERSION OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    list(APPEND Boost_b2_Command toolset=clang)
  elseif(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND Boost_b2_Command toolset=gcc)
  endif()
	
  if(ENV{CC})
    # CMake apprarently puts the full path of the compiler into CC
    # The user might specify a non-default gcc compiler through ENV
	message(STATUS "ENV{CC}=$ENV{CC}")
	get_filename_component( gccToolset "$ENV{CC}" NAME )

	# see: https://svn.boost.org/trac/boost/ticket/5917
	string(TOLOWER ${gccToolset} gccToolset)
	if(gccToolset STREQUAL "cc")
	  set(gccToolset "gcc")
	endif()
	list(APPEND Boost_b2_Command toolset=${gccToolset})
  endif()
  
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
	BUILD_COMMAND ${Boost_b2_Command} install -j8 --prefix=${Boost_Install_Dir} --with-thread --with-filesystem --with-system --with-date_time --with-program_options  --with-atomic  address-model=${Boost_address_model} link=static
	INSTALL_COMMAND ""
	)

  if(NOT WIN32)
    set(BOOST_ROOT        ${Boost_Install_Dir})
    set(BOOST_INCLUDE_DIR ${Boost_Install_Dir}/include)
  else()
    set(BOOST_ROOT        ${Boost_Install_Dir})
    set(Boost_INCLUDE_DIR ${Boost_Install_Dir}/include/boost-1_70)
  endif()

  set(Boost_LIBRARY_DIR ${Boost_Install_Dir}/lib)

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
