/*
 * feedbackMix.h
 *
 *  Created on: Aug 27, 2014
 *      Author: zhlou
 */

#ifndef FEEDBACKMIX_H_
#define FEEDBACKMIX_H_

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"
#include "aState.h"

template <class Problem>
class feedbackMix {
public:
    class Param {
    public:
        int interval;
        double target;
        Param(xmlNode *root);
    };
    feedbackMix(Problem &problem, const MPIState &mpiState,
                unirandom& rand, const Param &param);
    ~feedbackMix();
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
private:
    Mixing<Problem> mix;
    //xmlNode *root;
    double target;
    int interval;
    int tau_count;
    const MPIState &mpi;
    int *adoptArray;
    dynDebug debugOut;
};

#include "feedbackMix.hpp"

#endif /* FEEDBACKMIX_H_ */
