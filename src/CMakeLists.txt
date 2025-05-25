set(SFML_VERSION "3.0.0")
set(SFML_BUILD_EXAMPLES OFF CACHE BOOL "Disable SFML examples")

include(FetchContent)

set(SFML_VERSION 2.6.x)
FetchContent_Declare(
  SFML
  URL "https://github.com/SFML/SFML/archive/${SFML_VERSION}.zip"
)

## finally, let's build
FetchContent_MakeAvailable(sfml)



set(SOURCE_FILES
    main.cpp
    geometry_func.cpp
)
set(HEADER_FILES
    quad_tree.hpp
    geometry_func.hpp
    circle.hpp
    AABB.hpp
    particle.hpp
    particle.cpp
    free_list.hpp
)
add_executable(sfml-example
    ${SOURCE_FILES} ${HEADER_FILES}
)

target_link_libraries(sfml-example sfml-graphics )
target_include_directories(sfml-example PRIVATE "${SFML_INCLUDE_DIR}")
