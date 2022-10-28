add_executable(vector_indexing
  main.cpp
  )

target_link_libraries(vector_indexing
  Houdini
  ${GCOVLIB}
  )
