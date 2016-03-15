/*
 * exponential.cpp
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#include "exponential.h"
#include "xmlUtils.h"
#include <stdexcept>
#include <exception>

using namespace std;
const char *exponential::name = "exponential";


exponential::Param::Param(const ptree &root, debugStatus in_st, const char *name):
        st(in_st), outname(name)
{
    const ptree &section = root.get_child("exponential");
    alpha = section.get<double>("<xmlattr>.alpha");
    max_rej = section.get<unsigned long>("<xmlattr>.max_rej");
    segLength = section.get<unsigned>("<xmlattr>.log_freq",100);
}

exponential::exponential(const Param& param) : debugOut(param.st,param.outname), step_cnt(0)
{
    alpha = param.alpha;
    segLength = param.segLength;
    max_rej = param.alpha;
    reject_cnt = 0;
}


exponential::~exponential()
{
    // TODO Auto-generated destructor stub
}

void exponential::updateStep(bool accept, aState)
{
    if (accept)
        reject_cnt = 0;
    else
        reject_cnt ++;
}

void exponential::updateInitStep(bool, aState)
{

}

double exponential::updateS(aState state)
{
    return state.s / alpha;
}


void exponential::initStats(aState)
{
    reject_cnt = 0;
}


void exponential::updateSegment(aState state) {
    debugOut << state.step_cnt << "\t" << 1.0/state.s << "\t"
            << state.energy << std::endl;
}

void exponential::updateStats(aState state)
{
    step_cnt++;
    if (segLength == step_cnt) {
        step_cnt = 0;
        updateSegment(state);
        resetSegmentStats();
    }
}
