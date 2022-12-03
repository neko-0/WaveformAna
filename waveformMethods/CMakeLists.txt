cmake_minimum_required(VERSION 3.10)
project(WaveformMethods VERSION 1.0 DESCRIPTION "WaveformMethods")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(GNUInstallDirs)

add_library(WaveformMethods SHARED
  src/cfdTime.cpp
  src/baseline.cpp
  src/mix.cpp
  src/noise.cpp
  src/riseTime.cpp
  src/waveArea.cpp
  src/waveMax.cpp
  src/waveNMax.cpp
)

set_target_properties(WaveformMethods PROPERTIES
  VERSION ${PROJECT_VERSION}
  PUBLIC_HEADER waveformMethods/waveformMethods.hpp
)

target_include_directories(WaveformMethods PUBLIC .)

install(TARGETS WaveformMethods
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
