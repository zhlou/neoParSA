/*
 * expHold.cpp
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#include "expHold.h"
#include "xmlUtils.h"
#include <exception>
#include <stdexcept>

const char * expHold::name = "expHold";
expHold::Param::Param(xmlNode* root, debugStatus in_st, const char *name) :
        st(in_st), outname(name)
{
    xmlNode *xmlsection = getSectionByName(root, "expHold");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section expHold"));
    target_s = 1./getPropDouble(xmlsection, "target");
    alpha = getPropDouble(xmlsection, "alpha");
    segLength = 1;
    try {
        segLength = getPropInt(xmlsection, "segLength");
    } catch (const std::exception &e) {
        // ignored
    }
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
