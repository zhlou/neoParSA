/* 
 * File:   mixOnce.hpp
 * Author: zhlou
 *
 * Created on March 3, 2015, 3:52 PM
 */

#ifndef MIXONCE_HPP
#define	MIXONCE_HPP

#include "xmlUtils.h"

template <class Problem>
const char * mixOnce<Problem>::name = "mixOnce";

template <class Problem>
mixOnce<Problem>::Param::Param(xmlNode* root) : interval(1), useBest(true)
{
    xmlNode *section = getSectionByName(root, "mixOnce");
    target_s = 1./getPropDouble(section, "target");
    try {
        interval = getPropInt(section, "interval");
    } catch (const std::exception &e) {
        
    }
    try {
        useBest = getPropInt(section, "useBest");
    } catch (const std::exception &e)
    {
        
    }
    
}

template <class Problem>
mixOnce<Problem>::mixOnce(Problem& problem, const MPIState& mpiState, 
        unirandom& rnd, const Param& param) :
        mix(problem, mpiState, rnd),
        target_s(param.target_s), 
        interval(param.interval), 
        useBest(param.useBest),
        mixed(false),
        count(0)
{
    
}

template <class Problem>
mixState mixOnce<Problem>::Mix(aState& state)
{
    if (target_s > state.s)
        return mixState();
    if (mixed) {
        serialVar.update(state.energy);
        ++count;
        if ((count % interval) != 0)
            return mixState();
        count = 0;
        mix.calProbTab(state);
        debugOut << state.step_cnt << " " << mix.getEnergyVar() 
                << " " << serialVar.getVar() << '\n';
        return mixState();
    } else {
        mixed = true;
        count = 0;
        mix.calProbTab(state);
        int i;
        if (useBest)
            i = mix.getBest();
        else
            i = mix.getPartner();
        state.energy = mix.adoptState(i);
        return mixState(i);
        
    }
}

#endif	/* MIXONCE_HPP */

