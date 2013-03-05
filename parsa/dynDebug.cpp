/*
 * dynDebug.cpp
 *
 *  Created on: Mar 4, 2013
 *      Author: zhlou
 */
#include <iostream>
#include <stdexcept>
#include "dynDebug.h"

dynDebug::dynDebug(int i)
{
    if (i == 0) {
        isDebug = false;
        out = NULL;
    } else if (i == 1) {
        isDebug = true;
        out = &std::cout;
    } else if (i == 2) {
        isDebug = true;
        out = &std::cerr;
    } else {
        throw std::runtime_error("unrecognized debug option");
    }

}

dynDebug& dynDebug::operator<<(StandardEndLine)
{
    if (isDebug)
        *out << std::endl;
    return *this;
}
