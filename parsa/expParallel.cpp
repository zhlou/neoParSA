/*
 * pSimpleSchedule.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#include "expParallel.h"
#include "xmlUtils.h"
#include <exception>

const char *expParallel::name = "parallel exponential";

expParallel::expParallel(const Param &param, const MPIState &mpiState) :
    serExp(param.serParam), mpi(mpiState)
{
}

expParallel::Param::Param(xmlNode* root) : serParam(root)
{
}
