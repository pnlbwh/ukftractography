project(UKFCommonLib)

cmake_minimum_required(VERSION 2.8.7)

add_library(UKFCommonLib STATIC dummy.cc)

install(TARGETS UKFCommonLib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT Development
)

file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
install(FILES ${__files1} ${__files2}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include COMPONENT Development)

