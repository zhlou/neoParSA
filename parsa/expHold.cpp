/*
 * expHold.cpp
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#include "expHold.h"
#include <exception>
#include <stdexcept>

const char * expHold::name = "expHold";

expHold::Param::Param(const ptree &root, debugStatus in_st, const char *name) :
        st(in_st), outname(name)
{
    const ptree &section = root.get_child("expHold");
    target_s = 1./section.get<double>("<xmlattr>.target");
    alpha = section.get<double>("<xmlattr>.alpha");
    segLength = section.get<unsigned>("<xmlattr>.segLength", 1);
}

expHold::~expHold() {
    // TODO Auto-generated destructor stub
}

double expHold::updateS(const aState &state) {
    if (state.s < target_s)
        return state.s / alpha;
    else
        return target_s;
}

void expHold::updateSegment(const aState& state)
{
    debugOut << state.step_cnt << " " << state.s << " "
             << energyStat.getMean() << " "
             << energyStat.getVar() << std::endl;
}

void expHold::updateStats(const aState& state)
{
    step_cnt++;
    if (segLength == step_cnt) {
        step_cnt = 0;
        updateSegment(state);
        resetSegmentStats();
    }
}
