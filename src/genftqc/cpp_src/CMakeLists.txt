cmake_minimum_required(VERSION 3.14)
project(multiply_gates_libs)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

set(SOURCE_FILES
  cayley_graph.cpp
  cayley_graph.h)

set(WRAPPER_FILES
  cayley_graph_wrapper.cpp
  cayley_graph_wrapper.h
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

add_library(multiply_gates SHARED ${SOURCE_FILES} ${WRAPPER_FILES})
set_target_properties(multiply_gates PROPERTIES OUTPUT_NAME "multiply_gates")
target_link_libraries(multiply_gates ${HDF5_LIBRARIES})

if(WIN32)
  set_target_properties(multiply_gates PROPERTIES SUFFIX ".dll")
else()
  target_compile_options(multiply_gates PRIVATE -Wall -Wextra -pedantic)
endif()

add_executable(cayley_graph cayley_graph_main.cpp)
target_link_libraries(cayley_graph multiply_gates)
