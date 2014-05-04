set(TSP_COMMON tsp.cpp)
set(MPIBINARIES "")
set(BINARIES ${MPIBINARIES} solve_tsp test)
add_executable(solve_tsp solve_tsp.cpp ${TSP_COMMON})
add_executable(test test.cpp ${TSP_COMMON})
foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa ${LIBXML2_LIBRARIES} ${GSL_LIBRARIES})
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()


