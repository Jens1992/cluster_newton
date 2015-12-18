project(cluster_newton)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}")


set(COPASI_DEPENDENCY_DIR "" CACHE PATH "Dependency Directory containing, when set dependencies will be overwritten" )
mark_as_advanced(COPASI_DEPENDENCY_DIR)

if (COPASI_DEPENDENCY_DIR AND EXISTS ${COPASI_DEPENDENCY_DIR})
include(cmake-set-dependencies.cmake)
endif(COPASI_DEPENDENCY_DIR AND EXISTS ${COPASI_DEPENDENCY_DIR})


find_package(CLAPACK REQUIRED)

SET (SOURCES main.cpp)

include_directories(${CLAPACK_INCLUDE_DIR})
add_executable(cluster_newton ${SOURCES})
target_link_libraries(cluster_newton ${CLAPACK_LIBRARIES})
