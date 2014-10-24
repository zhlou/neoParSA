/*
 * adaptMix.hpp
 *
 *  Created on: Jan 31, 2013
 *      Author: zhlou
 */
#include <limits>
#include <cmath>
#include <iostream>
using namespace std;

template<class Problem>
const char *adaptMix<Problem>::name = "adaptiveMixing";

template<class Problem>
adaptMix<Problem>::Param::Param(xmlNode *root)
{
    xmlNode *section = getSectionByName(root, "mix");
    adaptCoef = getPropDouble(section, "adaptcoef");

}

template<class Problem>
adaptMix<Problem>::adaptMix(Problem& in_problem, const MPIState& mpiState,
        unirandom& in_rnd, const Param &param) :
        mix(in_problem, mpiState, in_rnd), rnd(in_rnd), nnodes(mpiState.nnodes),
        adaptCoef(param.adaptCoef)
{


}

template<class Problem>
adaptMix<Problem>::~adaptMix()
{

}

template<class Problem>
mixState adaptMix<Problem>::Mix(aState& state)
{
    mix.calProbTab(state);
    double rand = rnd.random();
    if (rand > adaptCoef * nnodes / mix.getNorm()) {
        int i = mix.getPartner();
        // cout << "rank " << mpi.rank << " adopts state from rank " << i << endl;
        state.energy = mix.adoptState(i);
        return (mixState(i)); // return (std::move(mixState(i)));
    } else {
        mix.adoptState();
        return mixState();
    }
}

