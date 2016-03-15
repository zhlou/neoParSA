set(NSCH_COMMON normSchwefel.cpp)
set(MPIBINARIES "")
set(BINARIES ${MPIBINARIES} nsch_scramble nsch_expHold)

add_executable(nsch_scramble nsch_scramble.cpp ${NSCH_COMMON})
add_executable(nsch_expHold nsch_expHold.cpp ${NSCH_COMMON})

foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa)
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()
add_custom_target(nsch DEPENDS ${BINARIES})
