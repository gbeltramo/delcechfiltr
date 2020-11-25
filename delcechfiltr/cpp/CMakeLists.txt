cmake_minimum_required(VERSION 3.11)
project(persty_cpp_binding VERSION 2.0.0
                           LANGUAGES CXX)

#============================
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
    set(PYBIND11_CPP_STANDARD -std=c++14)
    set(CMAKE_CXX_FLAGS "-O3")
endif()
if ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" )
    set(PYBIND11_CPP_STANDARD /std:c++14)
    set(CMAKE_CXX_FLAGS "/O2")
endif()

#============================
find_package(pybind11)

#============================
# UTIL LIBRARY
add_library(triangle STATIC
                     inc/triangle.hpp
                     src/triangle.cpp)
target_include_directories(triangle PUBLIC inc)
target_compile_features(triangle PUBLIC cxx_std_14)


#============================
# FINAL PYTHON MODULE USING PYBIND11
pybind11_add_module(binding MODULE binding.cpp)
target_link_libraries(binding PRIVATE triangle)
target_compile_features(binding PUBLIC cxx_std_14)
