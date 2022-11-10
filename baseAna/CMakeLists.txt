cmake_minimum_required(VERSION 3.10)

project(BaseAna VERSION 1.0 DESCRIPTION "BaseAna")

include(GNUInstallDirs)

add_library(BaseAna SHARED src/baseAna.cpp)

target_include_directories(BaseAna PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(BaseAna PUBLIC ConfigMgr)
