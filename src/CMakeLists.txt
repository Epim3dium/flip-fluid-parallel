set(SOURCE_FILES
    main.cpp
)
set(HEADER_FILES
)
add_executable(sfml-example
    ${SOURCE_FILES} ${HEADER_FILES}
)
find_package(SFML COMPONENTS graphics audio system REQUIRED)
target_link_libraries(sfml-example sfml-graphics sfml-window sfml-audio)
