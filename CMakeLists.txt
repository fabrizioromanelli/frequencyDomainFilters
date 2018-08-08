cmake_minimum_required(VERSION 2.8)

project(notch)
include_directories("include")
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11" )
set( CMAKE_VERBOSE_MAKEFILE ON )

set( DLS_FILTERS_LIB_NAME dlsfilters )

add_library(${DLS_FILTERS_LIB_NAME} SHARED src/bandFilter.cpp)

add_executable(${PROJECT_NAME} "src/notch.cpp")
target_link_libraries(${PROJECT_NAME} ${DLS_FILTERS_LIB_NAME})
