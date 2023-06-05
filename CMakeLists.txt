cmake_minimum_required(VERSION 3.16)

project(fv1d LANGUAGES C CXX)

add_subdirectory(external/HighFive)

include_directories(external/inih)

# Various sets
set(CMAKE_CXX_STANDARD 17)

add_executable(fv1d
  main.cpp)

target_link_libraries(fv1d PRIVATE HighFive)
