cmake_minimum_required(VERSION 3.10)

project(AnaQuickScan VERSION 1.0 DESCRIPTION "AnaQuickScan")

find_package(ROOT REQUIRED COMPONENTS RIO Net)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/AnaQuickScan)

include(GNUInstallDirs)

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIRS})

add_library(AnaQuickScan SHARED src/AnaQuickScan.cpp)

target_include_directories(AnaQuickScan PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaQuickScan PRIVATE BaseAna)
target_link_libraries(AnaQuickScan PRIVATE ConfigMgr WaveformMethods Utilities)
target_link_libraries(AnaQuickScan PUBLIC ROOT::Core ROOT::RIO ROOT::Net ${ROOT_LIBRARIES})
