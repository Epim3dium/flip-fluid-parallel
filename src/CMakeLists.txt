find_package(SFML 2.5.1 COMPONENTS system window graphics network audio REQUIRED)

set(USE_SEQ "seq")
set(USE_OPENMP "omp")
set(USE_CUDA "cuda")

set(BUILD_TYPE ${USE_CUDA})

set(SEQ_DIR "./sequential")
set(OPENMP_DIR "./openmp")
set(CUDA_DIR "./cuda")


set(SOURCE_FILES
    main.cpp
    geometry_func.cpp
)
set(HEADER_FILES
    geometry_func.hpp
    circle.hpp
    free_list.hpp
)
set(SEQ_SOURCE_FILES
    sequential/fluid.cpp
    sequential/particle.cpp
)
set(SEQ_HEADER_FILES
    sequential/fluid.hpp
    sequential/particle.hpp
    sequential/vec2.hpp
    sequential/AABB.hpp
)

set(OPENMP_SOURCE_FILES
    openmp/fluid.cpp
    openmp/particle.cpp
)
set(OPENMP_HEADER_FILES
    openmp/fluid.hpp
    openmp/particle.hpp
    openmp/vec2.hpp
    openmp/AABB.hpp
)
set(CUDA_SOURCE_FILES
    cuda/fluid.cpp
    cuda/particle.cu
    cuda/particle.cpp
)
set(CUDA_HEADER_FILES
    cuda/fluid.hpp
    cuda/particle.hpp
    cuda/vec2.hpp
    cuda/AABB.hpp
)
if("${BUILD_TYPE}" STREQUAL "${USE_SEQ}")
    message(STATUS "Using sequential code")
    add_executable(fluid-sim
        ${SOURCE_FILES} ${HEADER_FILES} ${SEQ_SOURCE_FILES} ${SEQ_HEADER_FILES}
    )
    target_include_directories(fluid-sim PRIVATE ${SEQ_DIR})
elseif("${BUILD_TYPE}" STREQUAL "${USE_OPENMP}")
    message(STATUS "Using openMP code")
    add_executable(fluid-sim
        ${SOURCE_FILES} ${HEADER_FILES} ${OPENMP_SOURCE_FILES} ${OPENMP_HEADER_FILES}
    )
    target_include_directories(fluid-sim PRIVATE ${OPENMP_DIR})
elseif("${BUILD_TYPE}" STREQUAL "${USE_CUDA}")
    message(STATUS "Using cuda code")
    add_executable(fluid-sim
        ${SOURCE_FILES} ${HEADER_FILES} ${CUDA_SOURCE_FILES} ${CUDA_HEADER_FILES}
    )
    target_include_directories(fluid-sim PRIVATE ${CUDA_DIR})
endif()
target_link_libraries(fluid-sim
  PUBLIC  sfml-graphics sfml-window sfml-system
)

target_link_libraries(fluid-sim PRIVATE sfml-graphics )
target_include_directories(fluid-sim PRIVATE "${SFML_INCLUDE_DIR}")
find_package(OpenMP REQUIRED)
target_link_libraries(fluid-sim PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(fluid-sim PRIVATE "${OpenMP_CXX_FLAGS}")
# target_compile_options(fluid-sim PRIVATE "${OpenMP_CXX_FLAGS}")

