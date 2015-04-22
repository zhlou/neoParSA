/*
 * intervalMix.hpp
 *
 *  Created on: Mar 21, 2013
 *      Author: zhlou
 */

#include "xmlUtils.h"
#include <limits>
#include <exception>
using namespace std;
//int readInterval(xmlNode *root);

template <class Problem>
const char * intervalMix<Problem>::name = "intervalMix";
/*
template<class Problem>
intervalMix<Problem>::intervalMix(Problem &in_problem, const MPIState &mpiState,
                                  unirandom& in_rand, xmlNode *docroot) :
        mix(in_problem, mpiState, in_rand), root(docroot),
        interval(readInterval(root)), tau_count(0)
{
}
*/

template<class Problem>
intervalMix<Problem>::intervalMix(Problem &in_problem, const MPIState &mpiState,
                                  unirandom& in_rand, const Param &param) :
        mix(in_problem, mpiState, in_rand),
        interval(param.interval), tau_count(0), reportNAdopt(param.reportNAdopt)
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
    if (reportNAdopt) {
        int nAdopt = mix.getNAdopt(i);
        debugOut << state.step_cnt << " " << state.s << " "
                << nAdopt << std::endl;
    }
    return mixState(i);
}

/*
inline int readInterval(xmlNode *root)
{
    xmlNode *section = getSectionByName(root, "mix");
    return getPropInt(section, "interval");
}
*/

template <class Problem>
intervalMix<Problem>::Param::Param(xmlNode *root)
{
    xmlNode *section = getSectionByName(root, "mix");
    interval = getPropInt(section, "interval");
    reportNAdopt = 0;
    try {
        reportNAdopt = getPropInt(section, "reportNAdopt");
    } catch (std::exception &e) {
        // ignore
    }
}

