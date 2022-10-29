add_executable(v2f_vector_indexing
  v2f_main.cpp
  )

target_link_libraries(v2f_vector_indexing
  Houdini
  ${GCOVLIB}
  )

Add_executable(struct_vector_indexing
  struct_main.cpp
  )

# gnuplot -e "plot 'plot_struct.txt' using 1:2 pt 7 ps 1 title 'Points', 'res_struct.txt' using 1:2 pt 7 ps 1 lc rgb 'red' title 'Result'; pause -1"
target_link_libraries(struct_vector_indexing
  Houdini
  ${GCOVLIB}
  )
