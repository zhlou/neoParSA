/* 
 * File:   rateMix.h
 * Author: zhlou
 *
 * Created on April 16, 2015, 1:02 PM
 */

#ifndef RATEMIX_H
#define	RATEMIX_H

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"
#include "aState.h"

template <class Problem>
class rateMix {
public:
    class Param {
    public:
        double factor;
        Param(xmlNode *root);
    };
    rateMix(Problem &problem, const MPIState &mpiState,
            unirandom &rand, const Param &param);
    ~rateMix();
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
private:
    Mixing<Problem> mix;
    const double factor;
    int interval;
    int tau_count;
    const MPIState &mpi;
    int *adoptArray;
    dynDebug debugOut;
};

#include "rateMix.hpp"

#endif	/* RATEMIX_H */

