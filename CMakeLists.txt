cmake_minimum_required(VERSION 2.8)

project(notch)
include_directories("include")
set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11" )
set( CMAKE_VERBOSE_MAKEFILE ON )

add_executable(${PROJECT_NAME} "src/notch.cpp")
