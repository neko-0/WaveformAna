cmake_minimum_required(VERSION 3.10)

project(WaveformAna VERSION 1.0 DESCRIPTION "WaveformAna")

add_subdirectory(configMgr)
add_subdirectory(utilities)
add_subdirectory(baseAna)
add_subdirectory(waveformMethods)
add_subdirectory(analysisDriver)

# user analysis directory
add_subdirectory(tests)
add_subdirectory(AnaSSRL)
