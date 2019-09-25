cmake_minimum_required(VERSION 2.8)

### -----------------------------------------------
### Project
### -----------------------------------------------
project(DG2D)

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

SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_executable(${PROJECT_NAME} ${ALLFILES})
include_directories(${dirlist})


### -----------------------------------------------
### MPI
### -----------------------------------------------

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
   set_target_properties(${PROJECT_NAME} PROPERTIES
   COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

### -----------------------------------------------
### OpenMP
### -----------------------------------------------
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


### -----------------------------------------------
### GTest
### -----------------------------------------------
add_subdirectory (ThirdParty/googletest)
enable_testing()

# Add test files
set( TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

add_executable( MeetGoogleTests ${TESTS_DIR}/simple_test.cpp )
target_link_libraries(MeetGoogleTests gtest gtest_main)
add_test( MeetGoogleTests MeetGoogleTests )

add_executable( CompServiceTests ${TESTS_DIR}/compServiceTests.cpp 
                                 ${CMAKE_CURRENT_SOURCE_DIR}/src/compService/compService.cpp
                                 ${CMAKE_CURRENT_SOURCE_DIR}/src/Mesh/Edge.cpp
                                 )
target_link_libraries(CompServiceTests gtest gtest_main)
add_test( CompServiceTests CompServiceTests )

### -----------------------------------------------
### END
### -----------------------------------------------
