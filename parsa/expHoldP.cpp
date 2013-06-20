/*
 * expHoldP.cpp
 *
 *  Created on: Jun 20, 2013
 *      Author: zhlou
 */

#include "expHoldP.h"

const char *expHoldP::name = "parallel expHold";
expHoldP::expHoldP(Param param, const MPIState &) : expHold(param.param)
{

}

expHoldP::~expHoldP() {
    // TODO Auto-generated destructor stub
}

