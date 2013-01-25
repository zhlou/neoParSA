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

simpleSchedule::simpleSchedule(xmlNode *root)
{
    xmlNode *xmlsection = getSectionByName(root, (char *)"annealer_input");

    if (xmlsection == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    //init_S = 1.0 / getPropDouble(xmlsection, (char *)"init_T");
    lambda = getPropDouble(xmlsection, (char *)"lambda");
    //init_loop = getPropInt(xmlsection, (char *)"init_loop");
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
    const unsigned max_rej = 100;
    return (reject_cnt >= max_rej);
}

double simpleSchedule::updateS(aState state)
{
    return state.s * (1 + lambda);
}

bool simpleSchedule::inSegment(aState)
{
    return false;
}

/*
double simpleSchedule::getInitS()
{
    return init_S;
}

int simpleSchedule::getInitLoop()
{
    return init_loop;
}
*/

void simpleSchedule::resetSegmentStats()
{
}

void simpleSchedule::updateSegment(aState)
{
}

void simpleSchedule::initStats(aState)
{
    reject_cnt = 0;
}
