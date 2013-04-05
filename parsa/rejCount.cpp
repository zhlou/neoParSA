/*
 * rejCount.cpp
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#include "rejCount.h"
#include "utils.h"

rejCount::rejCount(const Param &param) :
    debugOut(param.st, param.debugname)
{
    // TODO Auto-generated constructor stub
    reject_cnt = 0;
    max_rej = param.max_rej;
    // output_freq = param.output_freq;
}

rejCount::~rejCount() {
    // TODO Auto-generated destructor stub
}

void rejCount::updateState(bool accept, const aState &)
{
    if (accept)
        reject_cnt = 0;
    else
        ++ reject_cnt;
}

bool rejCount::frozen(const aState &state) const
{
    debugOut << state.step_cnt << "\t" << "reject_cnt" << std::endl;
    return (reject_cnt >= max_rej);
}

rejCount::Param rejCountParamXML(xmlNode *root, debugStatus st, const char *name)
{
    rejCount::Param param;
    param.st = st;
    param.debugname = name;
    xmlNode *xmlsection = getSectionByName(root, "count_reject");
    param.max_rej = getPropInt(xmlsection, "max_rej");

    return param;
}
