/*
 * tempCount.cpp
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#include "tempCount.h"

#include "xmlUtils.h"

tempCount::~tempCount() {
    // TODO Auto-generated destructor stub
}

tempCount::Param::Param(xmlNode* root, debugStatus in_st, const char* in_name)
{
    st = in_st;
    debugname = in_name;
    xmlNode *xmlsection = getSectionByName(root, "tempCount");
    target_s = 1.0/getPropDouble(xmlsection, "target");
    max_steps = getPropInt(xmlsection, "max_steps");
}

void tempCount::updateStep(bool, const aState& state) {
    if (state.s >= target_s)
        ++ step_cnt;
    else
        step_cnt = 0;
}

bool tempCount::frozen(const aState& state) {
    return ((state.s >= target_s) && (step_cnt >= max_steps));
}
