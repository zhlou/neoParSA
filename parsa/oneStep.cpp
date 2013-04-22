/*
 * oneStep.cpp
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#include "oneStep.h"
#include "utils.h"
#include <stdexcept>
#include <string>
const char * oneStep::name = "oneStep";
oneStep::oneStep(xmlNode *root) {
    xmlNode *xmlsection = getSectionByName(root, "oneStep");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section oneStep"));

    target_s = getPropDouble(xmlsection, "target");
}

oneStep::~oneStep() {
    // TODO Auto-generated destructor stub
}

