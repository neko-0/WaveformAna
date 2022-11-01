cmake_minimum_required(VERSION 3.10)

project(Utilities VERSION 1.0 DESCRIPTION "Utilities")

include(GNUInstallDirs)

add_library(Utilities SHARED
  src/logger.cpp
  src/getFiles.cpp
)

target_include_directories(Utilities PUBLIC ${CMAKE_CURRENT_LIST_DIR})
