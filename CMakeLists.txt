# MISANU -- koreni CMake: gradi OBA algoritma kao zasebne mete.
#   gvrp_pcafs_gvns  -> nas GVNS algoritam
#   mets_cpp         -> Xu-ov METS algoritam (C++ prevod)
# Otvori ovaj folder u CLion-u i dobijas obe mete u Run meniju.
cmake_minimum_required(VERSION 3.16)
project(MISANU CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_subdirectory(gvrp_pcafs_gvns)
add_subdirectory(mets_cpp)
