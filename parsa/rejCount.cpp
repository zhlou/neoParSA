/*
 * rejCount.cpp
 *
 *  Created on: Apr 3, 2013
 *      Author: zhlou
 */

#include "rejCount.h"

rejCount::rejCount(const Param &param) :
    debugOut(param.st, param.debugname)
{
    // TODO Auto-generated constructor stub
    reject_cnt = 0;
    max_rej = param.max_rej;
    output_freq = param.output_freq;
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
