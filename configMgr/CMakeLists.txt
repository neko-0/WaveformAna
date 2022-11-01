cmake_minimum_required(VERSION 3.10)

project(ConfigMgr VERSION 1.0 DESCRIPTION "ConfigMgr")

find_package(ROOT REQUIRED COMPONENTS RIO Net)

include(${ROOT_USE_FILE})
include(GNUInstallDirs)

# include_directories(./configMgr)

add_library(ConfigMgr SHARED
  src/betaConfigMgr.cpp
)
target_include_directories(ConfigMgr PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(ConfigMgr PUBLIC ROOT::RIO ROOT::Tree ROOT::Net)
target_link_libraries(ConfigMgr PUBLIC Utilities)

# set_target_properties(BetaConfigMgr PROPERTIES
#   VERSION ${PROJECT_VERSION}
#   PUBLIC_HEADER configMgr/betaConfigMgr.hpp
# )

# target_include_directories(BetaConfigMgr PUBLIC ./configMgr)

install(TARGETS ConfigMgr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
