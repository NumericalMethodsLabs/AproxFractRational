cmake_minimum_required(VERSION 3.17)
project(AproxFractRational)

set(CMAKE_CXX_STANDARD 17)

include_directories(MathMatrix)

add_subdirectory(MathMatrix)

add_library(${PROJECT_NAME} AproxFractRational.cpp AproxFractRational.h)

target_include_directories(${PROJECT_NAME} PUBLIC MathMatrix)

target_link_libraries(${PROJECT_NAME} MathMatrix)