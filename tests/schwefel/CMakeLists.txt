set(NSCH_COMMON normSchwefel.cpp)
set(MPIBINARIES "")
set(BINARIES ${MPIBINARIES} nsch_expHold)

add_executable(nsch_expHold nsch_expHold.cpp ${NSCH_COMMON})

foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa ${LIBXML2_LIBRARIES})
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()
add_custom_target(nsch DEPENDS ${BINARIES})
