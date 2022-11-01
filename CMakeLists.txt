cmake_minimum_required(VERSION 3.10)

project(WaveformAna VERSION 1.0 DESCRIPTION "WaveformAna")

add_subdirectory(configMgr)
add_subdirectory(analysisDriver)
add_subdirectory(waveformMethods)
add_subdirectory(utilities)

# user analysis directory
add_subdirectory(tests)
add_subdirectory(AnaSSRL)
