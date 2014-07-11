/*
 * expHoldP.h
 *
 *  Created on: Jun 20, 2013
 *      Author: zhlou
 */

#ifndef EXPHOLDP_H_
#define EXPHOLDP_H_

#include <mpi.h>
#include "MPIState.h"
#include "expHold.h"

class expHoldP: public expHold {
public:
    class Param{
    public:
        expHold::Param param;
        Param(xmlNode *root, debugStatus in_st=ignore, const char*name=NULL):
            param(root, in_st, name) {}
    };
    expHoldP(Param param, const MPIState &);
    ~expHoldP();
    static const char * name;
};

#endif /* EXPHOLDP_H_ */
