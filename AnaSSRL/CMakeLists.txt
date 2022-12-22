cmake_minimum_required(VERSION 3.10)

project(AnaSSRL VERSION 1.0 DESCRIPTION "AnaSSRL")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(AnaSSRL SHARED src/AnaSSRL.cpp)

target_include_directories(AnaSSRL PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_link_libraries(AnaSSRL PRIVATE BaseAna)
target_link_libraries(AnaSSRL PRIVATE ConfigMgr WaveformMethods Utilities)
