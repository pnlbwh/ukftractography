project(UKFCommonLib)

add_library(UKFCommonLib STATIC dummy.cc)

install(TARGETS UKFCommonLib
  RUNTIME DESTINATION bin COMPONENT RuntimeLibraries
  LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION lib COMPONENT Development
  )

file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
install(FILES ${__files1} ${__files2}
  DESTINATION include COMPONENT Development)

