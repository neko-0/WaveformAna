cmake_minimum_required(VERSION 3.10)

project(AnaQuickScan VERSION 1.0 DESCRIPTION "AnaQuickScan")

add_library(AnaQuickScan SHARED src/AnaQuickScan.cpp)

target_include_directories(AnaQuickScan PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaQuickScan PRIVATE BaseAna)
target_link_libraries(AnaQuickScan PRIVATE ConfigMgr WaveformMethods Utilities)
