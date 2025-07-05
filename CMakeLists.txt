cmake_minimum_required(VERSION 3.12)

set(USE_SEQ "seq")
set(USE_OPENMP "omp")
set(USE_CUDA "cuda")
set(BUILD_ALL "all")

set(BUILD_TYPE ${USE_OPENMP})

if("${BUILD_TYPE}" STREQUAL "${USE_CUDA}")
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
