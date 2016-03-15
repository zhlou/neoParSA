/*
 *       Filename:  staticCool.cpp
 *        Created:  09/21/2015 01:59:55 PM
 *         Author:  Zhihao Lou
 */

#include "staticCool.h"

#include <stdexcept>
#include <string>
#include "utils/vectorUtils.h"

const char *staticCool::name = "staticCool";


staticCool::Param::Param(const ptree &root, debugStatus in_st, const char*name):
    st(in_st), outname(name)
{
    const ptree &sec_attr = root.get_child("staticCool.<xmlattr>");
    std::string sname = sec_attr.get<std::string>("scheduleName");
    scheduleName = sname.c_str();
    segLength = sec_attr.get<unsigned>("segLength", 100);
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
