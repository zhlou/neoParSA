/*
 * periodBest.h
 *
 *  Created on: Aug 20, 2014
 *      Author: zhlou
 */

#ifndef PERIODBEST_H_
#define PERIODBEST_H_

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "unirandom.h"
#include "mixState.h"

template <class Problem>
class periodBest
{
private:
    Problem &problem;
    const MPIState &mpi;
    dynDebug debugOut;
    void *state_buf;
    const int buf_size;
    const unsigned frequency;
    unsigned counter;
    struct doubleint{
        double energy;
        int rank;
    };
public:
    class Param {
    public:
        unsigned frequency;
        Param(xmlNode *docroot);
    };
    periodBest(Problem &problem, MPIState const&mpiState,
               unirandom& in_rnd, const Param &param);
    ~periodBest();
    mixState Mix(aState &state);
    void setDebug(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st, outname);}
    static const char * name;
};

#include "periodBest.hpp"

#endif /* PERIODBEST_H_ */
