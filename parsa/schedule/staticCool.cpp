/*
 *       Filename:  staticCool.cpp
 *        Created:  09/21/2015 01:59:55 PM
 *         Author:  Zhihao Lou
 */

#include "staticCool.h"

#include <stdexcept>
#include <string>
#include "xmlUtils.h"
#include "utils/vectorUtils.h"

const char *staticCool::name = "staticCool";

staticCool::Param::Param(xmlNode *root, debugStatus in_st, const char*name):
    st(in_st), outname(name)
{
    xmlNode *xmlsection = getSectionByName(root, "staticCool");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section staticCool"));
    scheduleName = (char *)xmlGetProp(xmlsection, BAD_CAST"scheduleName");
    segLength=100;
    try {
        segLength = getPropInt(xmlsection, "segLength");
    } catch (const std::exception &e) {
        //ignored
    }
}

staticCool::staticCool(Param &param) :
    debugOut(param.st, param.outname),
    segLength(param.segLength),
    step_cnt(0),
    i(0)
{
    readVector(schedule, param.scheduleName);
    size = schedule.size();
}

staticCool::~staticCool()
{

}

double staticCool::updateS(const aState &state)
{
    if (i < size)
        return schedule[i++];
    else
        return schedule.back();
}


void staticCool::updateSegment(const aState &state)
{
    debugOut << state.step_cnt << " " << state.s << std::endl;
}

void staticCool::updateStats(const aState &state)
{
    ++ step_cnt;
    if (segLength == step_cnt) {
        step_cnt = 0;
        updateSegment(state);
        resetSegmentStats();
    }
}
