cmake_minimum_required(VERSION 3.12)

project(fluid-sim
  VERSION 1.0
  LANGUAGES CXX
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)


add_subdirectory(src)
