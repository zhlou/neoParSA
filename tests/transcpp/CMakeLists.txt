set(TRANSC_COMMON competition.cpp sequence.cpp score.cpp coeffects.cpp
                  cooperativity.cpp scalefactor.cpp fasta.cpp mode.cpp pwm.cpp
                  TF.cpp gene.cpp nuclei.cpp datatable.cpp
                  twobit.cpp parameter.cpp organism.cpp utils.cpp subgroup.cpp
                  bindings.cpp bindingsite.cpp distance.cpp promoter.cpp
                  quenching.cpp)
set(MPIBINARIES transcpp_pEF-mix
                transcpp_tempMove_pEF-mix
                transcpp_tempMove_pEF-fbMix
                transcpp_pEF-st-mix
                transcpp_pEF-fbMix
                transcpp_pEF-st-fbMix
                transcpp_pEF-twoTier
                transcpp_pEF-st-twoTier
                transcppPWLS)
#set(BINARIES ${MPIBINARIES} transc_lam transcpp_expHold transcpp_expFrozen
#             transc_scramble
#             transc_unfold parsefasta getsequence
#             transc_readFromState)

set(BINARIES transcpp transc_scramble transc_unfold parsefasta getsequence
             transc_readFromState
             transcpp_expFrozen
             transcpp_tempFB_expFrozen
             transcpp_expHold
             transcpp_lam
             transcpp_expHold-st
             transcpp_expFrozen-st
             ${MPIBINARIES})
add_executable(transcpp transcpp.cpp ${TRANSC_COMMON})
add_executable(transc_scramble scramble.cpp ${TRANSC_COMMON})
add_executable(transc_unfold unfold.cpp ${TRANSC_COMMON})
add_executable(parsefasta fasta.cpp parsefasta.cpp utils.cpp)
add_executable(getsequence twobit.cpp getsequence.cpp utils.cpp)

add_executable(transcpp_expHold transcpp_expHold.cpp ${TRANSC_COMMON})
add_executable(transcpp_expHold-st transcpp_expHold-st.cpp ${TRANSC_COMMON})
add_executable(transcpp_expFrozen transcpp_expFrozen.cpp ${TRANSC_COMMON})
add_executable(transcpp_tempFB_expFrozen transcpp_tempFB_expFrozen.cpp ${TRANSC_COMMON})
add_executable(transcpp_expFrozen-st transcpp_expFrozen-st.cpp ${TRANSC_COMMON})
add_executable(transcpp_lam transcpp_lam.cpp ${TRANSC_COMMON})
add_executable(transc_readFromState transc_readFromState.cpp ${TRANSC_COMMON})
#add_executable(transcpp_pEH-fbMix transcpp_pEH-fbMix.cpp ${TRANSC_COMMON})
#add_executable(transcpp_pEH-mix transcpp_pEH-mix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-mix transcpp_pEF-mix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-st-mix transcpp_pEF-st-mix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-fbMix transcpp_pEF-fbMix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-st-fbMix transcpp_pEF-st-fbMix.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-twoTier transcpp_pEF-twoTier.cpp ${TRANSC_COMMON})
add_executable(transcpp_pEF-st-twoTier transcpp_pEF-st-twoTier.cpp ${TRANSC_COMMON})
add_executable(transcppPWLS transcppPWLS.cpp ${TRANSC_COMMON})
add_executable(transcpp_tempMove_pEF-mix transcpp_tempMove_pEF-mix.cpp ${TRANSC_COMMON})
add_executable(transcpp_tempMove_pEF-fbMix transcpp_tempMove_pEF-fbMix.cpp ${TRANSC_COMMON})

foreach(bin ${BINARIES})
	target_link_libraries(${bin} parsa ${BOOST_LIBRARIES})
    #target_link_libraries(${bin} parsa ${BOOST_LIBRARIES})
endforeach()
foreach(mpibin ${MPIBINARIES})
	target_link_libraries(${mpibin} ${MPI_C_LIBRARIES})
endforeach()
add_custom_target(transc DEPENDS ${BINARIES})
