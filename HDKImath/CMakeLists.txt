cmake_minimum_required(VERSION 3.12)

project(HDKImath)

find_package(Houdini CONFIG REQUIRED)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS ( -pg -ftest-coverage -fprofile-arcs )
  SET (GCOVLIB "-lgcov")
endif()

add_executable(hdkimath
  main.cpp
  )

target_link_libraries(hdkimath
  Houdini
  ${GCOVLIB}
  )
