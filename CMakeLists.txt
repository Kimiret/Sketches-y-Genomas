cmake_minimum_required(VERSION 3.10)
project(Trabajo1)

add_subdirectory(smhasher)

add_executable(hyperloglog main.cpp)

target_link_libraries(hyperloglog PRIVATE smhasher)