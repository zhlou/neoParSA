/*
 * expHold.cpp
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#include "expHold.h"
#include "utils.h"
#include <stdexcept>

const char * expHold::name = "expHold";
expHold::expHold(xmlNode* root) {
    xmlNode *xmlsection = getSectionByName(root, "expHold");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section expHold"));
    target_s = 1./getPropDouble(xmlsection, "target");
    alpha = getPropDouble(xmlsection, "alpha");
}

expHold::~expHold() {
    // TODO Auto-generated destructor stub
}

double expHold::updateS(const aState &state) {
    if (state.s < target_s)
        return state.s / alpha;
    else
        return target_s;
}
