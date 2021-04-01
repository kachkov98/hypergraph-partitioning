cmake_minimum_required(VERSION 3.8)
project(partitioner)

set(SOURCES main.cpp FM.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

install (TARGETS ${PROJECT_NAME})
