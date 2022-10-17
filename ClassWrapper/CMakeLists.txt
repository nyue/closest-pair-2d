cmake_minimum_required(VERSION 3.12)

project(ClassWrapper)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS ( -pg -ftest-coverage -fprofile-arcs )
  SET (GCOVLIB "-lgcov")
endif()

add_executable(classwrapper
  main.cpp
  )

target_link_libraries(classwrapper ${GCOVLIB})
