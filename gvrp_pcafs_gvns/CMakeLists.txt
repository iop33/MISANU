cmake_minimum_required(VERSION 3.16)
project(gvrp_pcafs_gvns CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(gvrp_pcafs_gvns src/main.cpp)
target_include_directories(gvrp_pcafs_gvns PRIVATE src)
