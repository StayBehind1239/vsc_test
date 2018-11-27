cmake_minimum_required(VERSION 3.0.0)
project(vsc-test VERSION 0.1.0)

add_subdirectory(gsl)

include_directories(build/gsl)

add_executable(vsc-test main.cpp)

target_link_libraries(vsc-test gsl)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
