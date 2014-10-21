/*
 * dynDebug.cpp
 *
 *  Created on: Mar 4, 2013
 *      Author: zhlou
 */
#include <iostream>
#include <exception>
#include <stdexcept>
#include "dynDebug.h"

void dynDebug::__setDebug(const char* outname)
{
    if (status == ignore) {
        streamout = NULL;
    } else if (status == sout) {
        streamout = &std::cout;
    } else if (status == err) {
        streamout = &std::cerr;
    } else if (status == file) {
        try {
            streamout = new std::ofstream(outname);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "debug info is going to be ignored" << std::endl;
            status = ignore;
            streamout = NULL;
        }
    } else {
        throw std::runtime_error("unrecognized debug option");
    }
}

dynDebug::dynDebug(debugStatus st, const char *outname) : status(st)
{
    __setDebug(outname);
}

dynDebug& dynDebug::operator<<(StandardEndLine)
{
    if (status != ignore)
        *streamout << std::endl;
    return *this;
}

void dynDebug::setDebug(debugStatus st, const char *outname)
{
    if (status == file)
        delete streamout;
    status = st;
    __setDebug(outname);
}
