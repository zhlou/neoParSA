/* 
 * File:   mixOnce.h
 * Author: zhlou
 *
 * Created on March 3, 2015, 3:39 PM
 */

#ifndef MIXONCE_H
#define	MIXONCE_H

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"
#include "onePassMeanVar.h"

template <class Problem>
class mixOnce {
public:
    class Param {
    public:
        double target_s;
        int interval;
        bool useBest;
        Param(double target) : target_s(target), interval(1), useBest(true) {}
        Param(xmlNode *root);
    };
    mixOnce(Problem &problem, const MPIState &mpiState,
            unirandom& rnd, const Param &param);
    ~mixOnce(){}
    mixState Mix(aState &state);
    static const char * name;
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
private:
    Mixing<Problem> mix;
    const double target_s;
    const int interval;
    const bool useBest;
    bool mixed;
    int count;
    dynDebug debugOut;
    onePassMeanVar serialVar;
};

#include "mixOnce.hpp"

#endif	/* MIXONCE_H */

