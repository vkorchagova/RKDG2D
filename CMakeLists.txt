cmake_minimum_required(VERSION 2.8)

project(DG2D)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 ")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

FILE(GLOB ALLFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h
  )


MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${curdir}/${child})
    ENDIF()
  ENDFOREACH()
  LIST(APPEND dirlist ${curdir}/${child})
  SET(${result} ${dirlist})
ENDMACRO()


SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_executable(${PROJECT_NAME} ${ALLFILES})

include_directories(${dirlist})

target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
   set_target_properties(${PROJECT_NAME} PROPERTIES
   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
