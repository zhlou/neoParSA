/*
 * oneStep.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include "oneStep.h"
#include <stdexcept>
#include <string>
const char * oneStep::name = "oneStep";


oneStep::Param::Param(const ptree &root, debugStatus in_st, const char *name) :
        st(in_st), outname(name)
{
    target_s = 1./root.get<double>("oneStep.<xmlattr>.target");
}

oneStep::oneStep(Param param): target_s(param.target_s),
        debugOut(param.st, param.outname)
{}

oneStep::~oneStep() {
    // TODO Auto-generated destructor stub
}
