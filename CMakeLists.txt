cmake_minimum_required(VERSION 3.12)

option(USE_CUDA "Enable CUDA support" OFF)
option(USE_OMP "Enable OpenMP support" OFF)
option(BUILD_ALL "Build all" OFF)

if(USE_CUDA)
    project(fluid-sim
        VERSION 1.0
        LANGUAGES CXX C CUDA
    )
else()
    project(fluid-sim
        VERSION 1.0
        LANGUAGES CXX C
    )
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)


add_subdirectory(src)
