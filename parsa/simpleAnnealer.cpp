/*
 * simpleAnnealer.cpp
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#include "simpleAnnealer.h"

simpleAnnealer::simpleAnnealer(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    reject_cnt = 0;

}

simpleAnnealer::~simpleAnnealer()
{
    // TODO Auto-generated destructor stub
}

void simpleAnnealer::updateStep(bool accept, double delta)
{
    if (accept)
        reject_cnt = 0;
    else
        reject_cnt ++;
}

bool simpleAnnealer::frozen()
{
    const unsigned max_rej = 100;
    return (reject_cnt >= max_rej);
}

void simpleAnnealer::updateS()
{
    s *= (1 + lambda);
}

bool simpleAnnealer::inSegment()
{
    return false;
}

void simpleAnnealer::updateSegment()
{
}
