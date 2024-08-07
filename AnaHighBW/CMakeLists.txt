cmake_minimum_required(VERSION 3.10)

project(AnaHighBW VERSION 1.0 DESCRIPTION "AnaHighBW")

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(yaml-cpp REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIRS})

add_library(AnaHighBW SHARED src/AnaHighBW.cpp)

target_include_directories(AnaHighBW PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaHighBW PRIVATE BaseAna)
target_link_libraries(AnaHighBW PRIVATE ConfigMgr WaveformMethods Utilities)
target_link_libraries(AnaHighBW PRIVATE ${YAML_CPP_LIBRARIES})
