/*
 * intervalMix.hpp
 *
 *  Created on: Mar 21, 2013
 *      Author: zhlou
 */

#include "utils.h"
#include <limits>
using namespace std;
int readInterval(xmlNode *root);

template <class Problem>
const char * intervalMix<Problem>::name = "intervalMix";
template<class Problem>
intervalMix<Problem>::intervalMix(Problem &in_problem, const MPIState &mpiState,
                                  unirandom * const in_rand, xmlNode *docroot) :
        mix(in_problem, mpiState, in_rand), root(docroot),
        interval(readInterval(root)), tau_count(0)
{
}

template <class Problem>
intervalMix<Problem>::~intervalMix()
{
}

template <class Problem>
mixState intervalMix<Problem>::Mix(aState &state)
{
    ++tau_count;
    if ((tau_count % interval) != 0)
        return mixState();
    tau_count = 0;

    mix.calProbTab(state);
    int i = mix.getPartner();
    state.energy = mix.adoptState(i);
    return mixState(i);
}

inline int readInterval(xmlNode *root)
{
    xmlNode *section = getSectionByName(root, "mix");
    return getPropInt(section, "interval");
}
