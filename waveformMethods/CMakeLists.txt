cmake_minimum_required(VERSION 3.10)
project(WaveformMethods VERSION 1.0 DESCRIPTION "WaveformMethods")

include(GNUInstallDirs)

add_library(WaveformMethods SHARED
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
