/*
 * rejCount.cpp
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#include "rejCount.h"
#include "xmlUtils.h"

rejCount::rejCount(const rejCount::Param &param) :
    debugOut(param.st, param.debugname), max_rej(param.max_rej)
{
    // TODO Auto-generated constructor stub
    reject_cnt = 0;
    // output_freq = param.output_freq;
}


void rejCount::updateStep(bool accept, const aState &)
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

rejCount::Param::Param(xmlNode *root, debugStatus in_st, const char *name)
{
    st = in_st;
    debugname = name;
    xmlNode *xmlsection = getSectionByName(root, "count_reject");
    max_rej = getPropInt(xmlsection, "max_rej");
}
