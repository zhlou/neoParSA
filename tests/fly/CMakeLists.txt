find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

set(FLY_COMMON DataLists.cpp flyData.cpp maternal.cpp scoring.cpp solvers.cpp
	zygotic.cpp fly.cpp)
set(MPIBINARIES
    fly_sa-mpi
    fly_pEF-mix
    fly_pEF-fbMix
    fly_expHold-mpi
    fly_expHold-IndMove
    fly_expHold-theta
	fly_expHold-TA
    fly_expHold-mixing
    fly_expHold-mixing-IM
    fly_expHold-amix
	fly_simPLam
    fly_simPLam-fbMix
    fly_pEH-fbMix
    fly_pwls
    fly_pEF-tm-mix
    fly_pEF-tm-fbMix
    fly_cauchy_pEF-mix
    fly_cauchy_pEF-fbMix
)
set(BINARIES ${MPIBINARIES}
    printscore
    scramble
    fly_sa
    fly_pulse
    fly_expHold
    fly_expFrozen
    fly_expFrozen-tm
    fly_cauchy_expFrozen
    fly_wls
)

add_executable(printscore printscore.cpp ${FLY_COMMON})
add_executable(scramble scramble.cpp ${FLY_COMMON})
add_executable(fly_sa fly_sa.cpp ${FLY_COMMON})
add_executable(fly_pulse fly_pulse.cpp ${FLY_COMMON})
add_executable(fly_expHold fly_expHold.cpp ${FLY_COMMON})
add_executable(fly_expFrozen fly_expFrozen.cpp ${FLY_COMMON})
add_executable(fly_expFrozen-tm fly_expFrozen-tm.cpp ${FLY_COMMON})
add_executable(fly_sa-mpi fly_sa-mpi.cpp ${FLY_COMMON})
add_executable(fly_expHold-mpi fly_expHold-mpi.cpp ${FLY_COMMON})
add_executable(fly_expHold-IndMove fly_expHold-IndMove.cpp ${FLY_COMMON})
add_executable(fly_expHold-theta fly_expHold-theta.cpp ${FLY_COMMON})
add_executable(fly_expHold-TA fly_expHold-TA.cpp ${FLY_COMMON})
add_executable(fly_expHold-mixing fly_expHold-mixing.cpp ${FLY_COMMON})
add_executable(fly_expHold-mixing-IM fly_expHold-mixing-IM.cpp ${FLY_COMMON})
add_executable(fly_expHold-amix fly_expHold-amix.cpp ${FLY_COMMON})
add_executable(fly_simPLam fly_simPLam.cpp ${FLY_COMMON})
add_executable(fly_simPLam-fbMix fly_simPLam-fbMix.cpp ${FLY_COMMON})
add_executable(fly_pEH-fbMix fly_pEH-fbMix.cpp ${FLY_COMMON})
add_executable(fly_pEF-mix fly_pEF-mix.cpp ${FLY_COMMON})
add_executable(fly_pEF-fbMix fly_pEF-fbMix.cpp ${FLY_COMMON})
add_executable(fly_wls fly_wls.cpp ${FLY_COMMON})
add_executable(fly_pwls fly_pwls.cpp ${FLY_COMMON})
add_executable(fly_pEF-tm-mix fly_pEF-tm-mix.cpp ${FLY_COMMON})
add_executable(fly_pEF-tm-fbMix fly_pEF-tm-fbMix.cpp ${FLY_COMMON})
add_executable(fly_cauchy_expFrozen fly_cauchy_expFrozen.cpp ${FLY_COMMON})
add_executable(fly_cauchy_pEF-mix fly_cauchy_pEF-mix.cpp ${FLY_COMMON})
add_executable(fly_cauchy_pEF-fbMix fly_cauchy_pEF-fbMix.cpp ${FLY_COMMON})



foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa ${GSL_LIBRARIES})
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()
add_custom_target(fly DEPENDS ${BINARIES})
#set_target_properties(fly_sa-mpi PROPERTIES LINK_FLAGS "${MPI_C_LINK_FLAGS}")
