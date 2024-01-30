cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

project(WaveformAna VERSION 1.0 DESCRIPTION "WaveformAna")

add_subdirectory(configMgr)
add_subdirectory(utilities)
add_subdirectory(baseAna)
add_subdirectory(waveformMethods)
add_subdirectory(analysisDriver)

# user analysis directory
add_subdirectory(tests)
add_subdirectory(AnaSSRL)
add_subdirectory(AnaTCT)
add_subdirectory(AnaG4)
add_subdirectory(AnaSimple)
add_subdirectory(AnaQuickScan)
