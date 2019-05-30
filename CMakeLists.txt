cmake_minimum_required (VERSION 2.6)
project(cmdline_LBM)
include_directories("${CMAKE_SOURCE_DIR}/LBM" "${CMAKE_SOURCE_DIR}/LBM/headers" "${CMAKE_SOURCE_DIR}/rapidjson")
add_executable(main_cmdline ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp)
target_sources(main_cmdline
PRIVATE
  ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp
  ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp
  ${CMAKE_SOURCE_DIR}/rapidjson/document.h
  ${CMAKE_SOURCE_DIR}/rapidjson/writer.h
  ${CMAKE_SOURCE_DIR}/rapidjson/stringbuffer.h
)
target_link_libraries(main_cmdline)
