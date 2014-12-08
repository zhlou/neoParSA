find_package(Boost 1.52)
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DUSE_BOOST)
set(TRANSC_COMMON score.cpp coeffects.cpp cooperativity.cpp scalefactor.cpp 
                  fasta.cpp mode.cpp pwm.cpp TF.cpp gene.cpp nuclei.cpp conc.cpp 
                  twobit.cpp parameter.cpp organism.cpp utils.cpp subgroup.cpp 
                  bindings.cpp bindingsite.cpp distance.cpp promoter.cpp 
                  quenching.cpp)
set(MPIBINARIES transcpp_pEH-fbMix transcpp_pEH-mix)
set(BINARIES ${MPIBINARIES} transc_lam transcpp_expHold transc_scramble 
             transc_unfold parsefasta getsequence)

add_executable(transc_lam transcpp.cpp ${TRANSC_COMMON})
add_executable(transcpp_expHold transcpp_expHold.cpp ${TRANSC_COMMON})
add_executable(transc_scramble scramble.cpp ${TRANSC_COMMON})
add_executable(transc_unfold unfold.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEH-fbMix transcpp_pEH-fbMix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEH-mix transcpp_pEH-mix.cpp ${TRANSC_COMMON})
add_executable(parsefasta fasta.cpp parsefasta.cpp)
add_executable(getsequence twobit.cpp getsequence.cpp)

foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa ${LIBXML2_LIBRARIES} ${BOOST_LIBRARIES})
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()
add_custom_target(transcpp DEPENDS ${BINARIES})
