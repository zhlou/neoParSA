/*
 * tempCountP.h
 *
 *  Created on: Jun 20, 2013
 *      Author: zhlou
 */

#ifndef TEMPCOUNTP_H_
#define TEMPCOUNTP_H_

#include "tempCount.h"
#include "MPIState.h"
#include <mpi.h>

class tempCountP: public tempCount {
public:
    class Param {
    public:
        tempCount::Param param;
        Param(xmlNode *root, debugStatus st=ignore, const char * debugname = NULL) :
            param(root, st, debugname) {}
    };
    tempCountP(const Param &param, const MPIState &mpiState);
    virtual ~tempCountP();
    // static const char * name;
};

#endif /* TEMPCOUNTP_H_ */
