/*
 * twoTierMix.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: zhlou
 */

#ifndef TWOTIERMIX_HPP_
#define TWOTIERMIX_HPP_

#include "xmlUtils.h"

template <class Problem>
twoTierMix<Problem>::Param::Param(xmlNode *docroot)
{
    xmlNode *section = getSectionByName(docroot, "twoTierMix");
    partSize = getPropInt(section, "partSize");
    mixFreq = getPropInt(section, "mixFreq");
    globalFreq = getPropInt(section, "globalFreq");
}

template <class Problem>
const char *twoTierMix<Problem>::name = "twoTierMix";

template <class Problem>
twoTierMix<Problem>::twoTierMix(
        Problem &in_problem,
        const MPIState &globalState,
        unirandom& in_rand,
        const Param &param):
            globalState(globalState),
            globalMix(in_problem, globalState, in_rand),
            mixFreq(param.mixFreq),
            globalFreq(param.globalFreq),
            tauCount(0),
            mixCount(0)
{
    MPI_Comm_split(globalState.comm, globalState.rank / param.partSize,
                   param.partSize, &localState.comm);
    MPI_Comm_size(localState.comm, &localState.nnodes);
    MPI_Comm_rank(localState.comm, &localState.rank);
    MPI_Comm_group(localState.comm, &localState.group);
    localMix = new Mixing<Problem>(in_problem, localState, in_rand);
    local2Global = new int[localState.nnodes];
    int *seq = new int[localState.nnodes];
    for (int i = 0; i < localState.nnodes; ++i)
        seq[i] = i;
    MPI_Group_translate_ranks(localState.group, localState.nnodes,
                              seq, globalState.group, local2Global);
    delete []seq;


}

template <class Problem>
twoTierMix<Problem>::~twoTierMix()
{
    delete localMix;
    MPI_Group_free(&localState.group);
    MPI_Comm_free(&localState.comm);
    delete []local2Global;
}

template<class Problem>
mixState twoTierMix<Problem>::Mix(aState &state)
{
    ++tauCount;
    if ((tauCount % mixFreq) != 0)
        return mixState();
    tauCount = 0;
    ++mixCount;
    if ((mixCount % globalFreq) != 0) {
        localMix->calProbTab(state);
        int i = localMix->getPartner();
        state.energy = localMix->adoptState(i);
        return mixState(local2Global[i]);
    } else {
        mixCount = 0;
        globalMix.calProbTab(state);
        int i = globalMix.getPartner();
        state.energy = globalMix.adoptState(i);
        return mixState(i);
    }
}

#endif /* TWOTIERMIX_HPP_ */
