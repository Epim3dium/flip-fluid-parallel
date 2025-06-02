find_package(SFML 2.5.1 COMPONENTS system window graphics network audio REQUIRED)


option(USE_SEQ "Use sequential backend" OFF)
option(USE_OPENMP "Use OpenMP backend" ON)
option(USE_CUDA "Use CUDA backend" OFF)
set(SEQ_DIR "./sequential")
set(OPENMP_DIR "./openmp")

set(SOURCE_FILES
    main.cpp
    geometry_func.cpp
)
set(HEADER_FILES
    geometry_func.hpp
    circle.hpp
    AABB.hpp
    free_list.hpp
)
set(SEQ_SOURCE_FILES
    sequential/fluid.cpp
    sequential/particle.cpp
)
set(SEQ_HEADER_FILES
    sequential/fluid.hpp
    sequential/particle.hpp
)

set(OPENMP_SOURCE_FILES
    openmp/fluid.cpp
    openmp/particle.cpp
)
set(OPENMP_HEADER_FILES
    openmp/fluid.hpp
    openmp/particle.hpp
)
if(${USE_SEQ})
    add_executable(fluid-sim
        ${SOURCE_FILES} ${HEADER_FILES} ${SEQ_SOURCE_FILES} ${SEQ_HEADER_FILES}
    )
    target_include_directories(fluid-sim PRIVATE ${SEQ_DIR})
elseif(${USE_OPENMP})
    add_executable(fluid-sim
        ${SOURCE_FILES} ${HEADER_FILES} ${OPENMP_SOURCE_FILES} ${OPENMP_HEADER_FILES}
    )
    target_include_directories(fluid-sim PRIVATE ${OPENMP_DIR})
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

