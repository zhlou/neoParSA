/*
 * simpleAnnealer.cpp
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#include "simpleScheduler.h"
#include "utils.h"
#include <stdexcept>

using namespace std;
const char *simpleSchedule::name = "exponential";

simpleSchedule::simpleSchedule(xmlNode *root)
{
    xmlNode *xmlsection = getSectionByName(root, "exponential");

    if (xmlsection == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    //init_S = 1.0 / getPropDouble(xmlsection, "init_T");
    alpha = getPropDouble(xmlsection, "alpha");
    max_rej = getPropInt(xmlsection, "max_rej");
    //init_loop = getPropInt(xmlsection, "init_loop");
    reject_cnt = 0;

}

simpleSchedule::~simpleSchedule()
{
    // TODO Auto-generated destructor stub
}

void simpleSchedule::updateStep(bool accept, aState)
{
    if (accept)
        reject_cnt = 0;
    else
        reject_cnt ++;
}

void simpleSchedule::updateInitStep(bool, aState)
{

}

bool simpleSchedule::frozen(aState)
{
    // const unsigned max_rej = 100;
    return (reject_cnt >= max_rej);
}

double simpleSchedule::updateS(aState state)
{
    return state.s * (1 + alpha);
}


void simpleSchedule::initStats(aState)
{
    reject_cnt = 0;
}

void simpleSchedule::updateSegment(aState state) {
    debugOut << state.step_cnt << "\t" << 1.0/state.s << "\t"
            << state.energy << std::endl;
}
