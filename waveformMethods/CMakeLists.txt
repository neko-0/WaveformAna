cmake_minimum_required(VERSION 3.10)
project(WaveformMethods VERSION 1.0 DESCRIPTION "WaveformMethods")

find_package(OpenMP)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include(GNUInstallDirs)

add_library(WaveformMethods SHARED
  src/baseline.cpp
  src/cfdTime.cpp
  src/mix.cpp
  src/noise.cpp
  src/riseTime.cpp
  src/utils.cpp
  src/waveArea.cpp
  src/waveMax.cpp
  src/waveNMax.cpp
)

set(WaveformMethods_HEADERS
  "waveformMethods/core.hpp"
  "waveformMethods/baseline.hpp"
  "waveformMethods/utils.hpp"
  "waveformMethods/waveformMethods.hpp"
)

set_target_properties(WaveformMethods PROPERTIES
  VERSION ${PROJECT_VERSION}
  PUBLIC_HEADER "${WaveformMethods_HEADERS}"
)

target_include_directories(WaveformMethods PUBLIC .)

target_link_libraries(WaveformMethods Eigen3::Eigen)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    target_link_libraries (WaveformMethods OpenMP::OpenMP_CXX)
endif()

install(TARGETS WaveformMethods
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
