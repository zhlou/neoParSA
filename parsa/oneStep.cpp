/*
 * oneStep.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include "oneStep.h"
#include "xmlUtils.h"
#include <stdexcept>
#include <string>
const char * oneStep::name = "oneStep";

oneStep::Param::Param(xmlNode *root, debugStatus in_st, const char *name) :
        st(in_st), outname(name)
{
    xmlNode *xmlsection = getSectionByName(root, "oneStep");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section oneStep"));

    target_s = 1./getPropDouble(xmlsection, "target");
}

oneStep::oneStep(Param param): target_s(param.target_s),
        debugOut(param.st, param.outname)
{}

oneStep::~oneStep() {
    // TODO Auto-generated destructor stub
}

