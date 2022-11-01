cmake_minimum_required(VERSION 3.10)

project(Tester VERSION 1.0 DESCRIPTION "Tester")

add_executable(tester tester.cpp)
# set_target_properties(EXEC PROPERTIES OUTPUT_NAME "EXEC")
target_link_libraries(tester PRIVATE ConfigMgr WaveformMethods)
