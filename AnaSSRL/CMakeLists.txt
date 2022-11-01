cmake_minimum_required(VERSION 3.10)

project(AnaSSRL VERSION 1.0 DESCRIPTION "AnaSSRL")

add_library(AnaSSRL SHARED src/AnaSSRL.cpp)

target_include_directories(AnaSSRL PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaSSRL PRIVATE AnalysisDriver)
target_link_libraries(AnaSSRL PRIVATE ConfigMgr WaveformMethods Utilities)
