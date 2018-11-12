cmake_minimum_required(VERSION 2.8)

project(decomposeUNV)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/*.h
  )

add_executable(${PROJECT_NAME} ${ALLFILES})
