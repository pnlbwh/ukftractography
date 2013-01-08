# -----------------------------------------------------------------------------
# GETSLICERVERSION
#   - Tries to find the different Slicer version starting with the latest
#   - When Slicer is found it includes it
#   - Slicer_VERSION is an output variable
#-----------------------------------------------------------------------------

if(NOT GETSLICERVERSION)
macro(GETSLICERVERSION Slicer_VERSION)
  
  # Slicer4 complains if those variables are set during ccmake
  unset(ITK_DIR CACHE)
  unset(CTK_DIR CACHE)
  unset(Teem_DIR CACHE)

  # Find Slicer4
  find_package(Slicer QUIET)

  if (Slicer_FOUND)

    include(${Slicer_USE_FILE})
    set(Slicer_VERSION 4)

  else (Slicer_FOUND)

    # if not found find Slicer3
    find_package(Slicer3 QUIET)
   
  endif (Slicer_FOUND)

  if (Slicer3_FOUND)
    include(${Slicer3_USE_FILE})
    set(Slicer_VERSION 3)
  endif (Slicer3_FOUND)

  # If no Slicer version is found set the version to zero
  if (NOT Slicer_FOUND AND NOT Slicer3_FOUND)
    set(Slicer_VERSION 0)
  endif (NOT Slicer_FOUND AND NOT Slicer3_FOUND)

endmacro(GETSLICERVERSION Slicer_VERSION)
endif(NOT GETSLICERVERSION)

#-----------------------------------------------------------------------------
# PARSE_ARGUMENTS
#   - Parses the arguments of a function
#   - See www.cmake.org/Wiki/CMakeMacroParseArguments for details
#-----------------------------------------------------------------------------

IF( NOT PARSE_ARGUMENTS)
MACRO(PARSE_ARGUMENTS prefix arg_names option_names)

 SET(DEFAULT_ARGS)
 FOREACH(arg_name ${arg_names})   
   SET(${prefix}_${arg_name})
 ENDFOREACH(arg_name)
 FOREACH(option ${option_names})
   SET(${prefix}_${option} FALSE)
 ENDFOREACH(option)

  SET(current_arg_name DEFAULT_ARGS)
  SET(current_arg_list)
  FOREACH(arg ${ARGN})           
    SET(larg_names ${arg_names})   
    LIST(FIND larg_names "${arg}" is_arg_name)                   
    IF (is_arg_name GREATER -1)
      SET(${prefix}_${current_arg_name} ${current_arg_list})
      SET(current_arg_name ${arg})
      SET(current_arg_list)
    ELSE (is_arg_name GREATER -1)
      SET(loption_names ${option_names})   
      LIST(FIND loption_names "${arg}" is_option)           
      IF (is_option GREATER -1)
             SET(${prefix}_${arg} TRUE)
      ELSE (is_option GREATER -1)
             SET(current_arg_list ${current_arg_list} ${arg})
      ENDIF (is_option GREATER -1)
    ENDIF (is_arg_name GREATER -1)
  ENDFOREACH(arg)
  SET(${prefix}_${current_arg_name} ${current_arg_list})

ENDMACRO(PARSE_ARGUMENTS) 
ENDIF(NOT PARSE_ARGUMENTS)

#-----------------------------------------------------------------------------
# GENERATECLPMACRO
#   - Runs GENERATECLP with the appropriate variables
#   - Mimics slicerMacroBuildCLI of Slicer4
#-----------------------------------------------------------------------------
if(NOT GENERATECLPMACRO)
macro(GENERATECLPMACRO)

  PARSE_ARGUMENTS(ARG "ADDITIONAL_SRCS;XML_FILE;EXECUTABLE_NAME;TARGET_LIBRARIES;INCLUDE_DIRECTORIES;MAIN_SRC" "" ${ARGN})

  ## Additional include paths, link required for boost 
  if (${ARG_INCLUDE_DIRECTORIES}) 
    include(${ARG_INCLUDE_DIRECTORIES})
  endif (${ARG_INCLUDE_DIRECTORIES}) 

  ## Generate InterfaceCLP.h
  ## Macro defined in SlicerExecutionModel
  GENERATECLP(${ARG_MAIN_SRC} ${ARG_XML_FILE})

  ## Add Executables and Link
  ADD_EXECUTABLE(${ARG_EXECUTABLE_NAME} ${ARG_MAIN_SRC} ${ARG_ADDITIONAL_SRCS} )
  TARGET_LINK_LIBRARIES(${ARG_EXECUTABLE_NAME} ${ARG_TARGET_LIBRARIES})

endmacro(GENERATECLPMACRO)
endif(NOT GENERATECLPMACRO)

#-----------------------------------------------------------------------------
# SETIFEMPTY
#   - Sets a variable only if it's not already set
#-----------------------------------------------------------------------------
if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${ARGV})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)
