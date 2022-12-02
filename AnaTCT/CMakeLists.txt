cmake_minimum_required(VERSION 3.10)

project(AnaTCT VERSION 1.0 DESCRIPTION "AnaTCT")

add_library(AnaTCT SHARED src/AnaTCT.cpp)

target_include_directories(AnaTCT PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaTCT PRIVATE BaseAna)
target_link_libraries(AnaTCT PRIVATE ConfigMgr WaveformMethods Utilities)
