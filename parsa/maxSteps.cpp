/*
 * maxSteps.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include "maxSteps.h"
#include "xmlUtils.h"
maxSteps::Param::Param(xmlNode *root, debugStatus in_st, const char *name)
{
    st = in_st;
    debugname = name;
    xmlNode *xmlsection = getSectionByName(root, "maxSteps");
    max_steps = getPropInt(xmlsection, "max_steps");
}

maxSteps::~maxSteps() {
    // TODO Auto-generated destructor stub
}

