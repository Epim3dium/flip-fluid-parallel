find_package(SFML 2.5.1 COMPONENTS system window graphics network audio REQUIRED)



set(SOURCE_FILES
    main.cpp
    geometry_func.cpp
    fluid.cpp
)
set(HEADER_FILES
    geometry_func.hpp
    circle.hpp
    AABB.hpp
    particle.hpp
    particle.cpp
    free_list.hpp
    fluid.hpp
)
add_executable(fluid-sim
    ${SOURCE_FILES} ${HEADER_FILES}
)
target_link_libraries(fluid-sim
  PUBLIC  sfml-graphics sfml-window sfml-system
)

target_link_libraries(fluid-sim PRIVATE sfml-graphics )
target_include_directories(fluid-sim PRIVATE "${SFML_INCLUDE_DIR}")
find_package(OpenMP REQUIRED)
target_link_libraries(fluid-sim PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(fluid-sim PRIVATE "${OpenMP_CXX_FLAGS}")
# target_compile_options(fluid-sim PRIVATE "${OpenMP_CXX_FLAGS}")

