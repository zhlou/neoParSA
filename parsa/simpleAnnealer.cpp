/*
 * simpleAnnealer.cpp
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#include "simpleAnnealer.h"
#include <exception>

using namespace std;

simpleSchedule::simpleSchedule(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    xmlNode *xmlsection = getSectionByName(root, "annealer_input");

    if (xmlsection == NULL) {
        throw runtime_error(string("Error: fail to find section annealer_input"));
    }
    init_S = 1.0 / getPropDouble(xmlsection, "init_T");
    lambda = getPropDouble(xmlsection, "lambda");
    init_loop = getPropInt(xmlsection, "init_loop");
    reject_cnt = 0;

}

simpleSchedule::~simpleSchedule()
{
    // TODO Auto-generated destructor stub
}

void simpleSchedule::updateStep(bool accept, double)
{
    if (accept)
        reject_cnt = 0;
    else
        reject_cnt ++;
}

inline void simpleSchedule::updateInitStep(bool, double)
{

}

bool simpleSchedule::frozen()
{
    const unsigned max_rej = 100;
    return (reject_cnt >= max_rej);
}

double simpleSchedule::updateS(double s)
{
    return s * (1 + lambda);
}

bool simpleSchedule::inSegment()
{
    return false;
}

inline double simpleSchedule::getInitS()
{
    return init_S;
}

inline int simpleSchedule::getInitLoop()
{
    return init_loop;
}

void simpleSchedule::resetSegmentStats()
{
}

void simpleSchedule::updateSegment()
{
}

void simpleSchedule::initStats()
{
    reject_cnt = 0;
}
