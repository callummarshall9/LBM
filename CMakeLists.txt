cmake_minimum_required (VERSION 2.6)
include_directories("LBM")
add_subdirectory(LBM)
include_directories("LBM/headers")
add_subdirectory(LBM/headers)
project (cmdline_LBM)
add_executable(main_cmdline LBM/source/LBM.cpp LBM/source/vector3.cpp cmdline/main.cpp)

project(gui_LBM)
find_package(PkgConfig)
find_package (Threads)
pkg_check_modules(GTKMM gtkmm-3.0)

link_directories(
    ${GTKMM_LIBRARY_DIRS}  )

include_directories(
  ${GTKMM_INCLUDE_DIRS} "LBM" "LBM/headers" "gui/headers" )
add_executable(main_gui)
target_sources(main_gui
PRIVATE
  gui/main.cpp
  LBM/source/LBM.cpp
  LBM/source/vector3.cpp
  gui/source/worker.cpp
)
target_link_libraries(main_gui Threads::Threads ${GTKMM_LIBRARIES})
