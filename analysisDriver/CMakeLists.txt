cmake_minimum_required(VERSION 3.10)

project(AnalysisDriver VERSION 1.0 DESCRIPTION "AnalysisDriver")

# find_package(ROOT REQUIRED COMPONENTS RIO Net)
# include(${ROOT_USE_FILE})

include(GNUInstallDirs)

add_library(AnalysisDriver SHARED
  src/register.cpp
  src/driver.cpp
)

# Need to include user analysis here as well.
set(ANA_LIST AnaSSRL AnaTCT)

target_include_directories(AnalysisDriver PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_link_libraries(AnalysisDriver PUBLIC BaseAna)
target_link_libraries(AnalysisDriver PUBLIC ConfigMgr)
target_link_libraries(AnalysisDriver PUBLIC Utilities)
target_link_libraries(AnalysisDriver PUBLIC ${ANA_LIST})
# target_link_libraries(AnalysisDriver PUBLIC ROOT::RIO ROOT::Tree ROOT::Net)

# install(TARGETS AnalysisDriver
#     LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#     PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
# )

#===============================================================================
project(RunAna VERSION 1.0 DESCRIPTION "RunAna")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(CMAKE_CXX_FLAGS "-lboost_program_options")
# find_package(Boost REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(run_Ana utils/run_Ana.cpp)
    target_link_libraries(run_Ana PUBLIC ${Boost_LIBRARIES})
    target_link_libraries(run_Ana PUBLIC AnalysisDriver ConfigMgr Utilities)
endif()

add_executable(ls_Ana utils/ls_Ana.cpp)
target_link_libraries(ls_Ana PUBLIC AnalysisDriver Utilities)

add_executable(check_Ana utils/check_Ana.cpp)
target_link_libraries(check_Ana PUBLIC AnalysisDriver Utilities)
