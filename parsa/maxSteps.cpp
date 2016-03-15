/*
 * maxSteps.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include "maxSteps.h"


maxSteps::Param::Param(const ptree &root, debugStatus in_st, const char *name):
        st(in_st), debugname(name)
{
    max_steps = root.get<int>("maxSteps.<xmlattr>.max_steps");
}

maxSteps::~maxSteps() {
    // TODO Auto-generated destructor stub
}
