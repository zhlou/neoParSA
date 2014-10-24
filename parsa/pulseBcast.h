/*
 * pulseBcast.h
 *
 *  Created on: May 15, 2013
 *      Author: zhlou
 */

#ifndef PULSEBCAST_H_
#define PULSEBCAST_H_

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "unirandom.h"
#include "mixState.h"

template <class Problem>
class pulseBcast
{
private:
    Problem &problem;
    const MPIState &mpi;
    unirandom& rnd;
    void *state_buf;
    const int buf_size;
    unsigned counter;
    const unsigned frequency;
    dynDebug debugOut;
public:
    class Param {
    public:
        unsigned frequency;
        Param(xmlNode *root);
    };
    pulseBcast(Problem &problem, const MPIState &mpiState,
               unirandom& in_rnd, const Param &param);
    ~pulseBcast();
    mixState Mix(aState &state);
    void setDebug(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st,outname);}
    static const char * name;

};

#include "pulseBcast.hpp"

#endif /* PULSEBCAST_H_ */
