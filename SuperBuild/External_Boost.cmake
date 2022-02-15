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
    set(Boost_url "http://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.tar.gz")
    set(Boost_SHA256 94ced8b72956591c4775ae2207a9763d3600b30d9d7446562c552f0a14a63be7)
    set(Boost_Bootstrap_Command ./bootstrap.sh)
    set(Boost_b2_Command ./b2)
  else()
    if(WIN32)
      set(Boost_url "http://sourceforge.net/projects/boost/files/boost/1.78.0/boost_1_78_0.zip")
      set(Boost_SHA256 f22143b5528e081123c3c5ed437e92f648fe69748e95fa6e2bd41484e2986cc3)
      set(Boost_Bootstrap_Command ./bootstrap.bat)
      set(Boost_b2_Command b2)
    endif()
  endif()

  if(MSVC)
    if(MSVC_VERSION GREATER_EQUAL 1800 AND MSVC_VERSION LESS 1900)
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
  unset(deplyment_mode)
  if(APPPLE)
     set(deployment_mode "-a macosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
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
	URL_HASH SHA256=${Boost_SHA256}
	UPDATE_COMMAND ""
	CONFIGURE_COMMAND ${Boost_Bootstrap_Command} --prefix=${Boost_Install_Dir}/lib
	BUILD_COMMAND ${Boost_b2_Command} ${deployment_mode} install -j8 --prefix=${Boost_Install_Dir} --with-thread --with-filesystem --with-system --with-date_time --with-program_options  --with-atomic  address-model=${Boost_address_model} link=static
	INSTALL_COMMAND ""
	)

  if(NOT WIN32)
    set(BOOST_ROOT        ${Boost_Install_Dir})
    set(BOOST_INCLUDE_DIR ${Boost_Install_Dir}/include)
  else()
    set(BOOST_ROOT        ${Boost_Install_Dir})
    set(Boost_INCLUDE_DIR ${Boost_Install_Dir}/include/boost-1_78)
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
