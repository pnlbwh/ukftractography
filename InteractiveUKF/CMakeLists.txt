
#-----------------------------------------------------------------------------
set(MODULE_NAME InteractiveUKF)
set(MODULE_TITLE ${MODULE_NAME})

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

#-----------------------------------------------------------------------------
add_subdirectory(Logic)
add_subdirectory(UI)

#-----------------------------------------------------------------------------
set(MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT")

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  #  add_subdirectory(Testing)
endif()
