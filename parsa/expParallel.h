/*
 * pSimpleSchedule.h
 *
 *  Created on: Apr 2, 2013
 *      Author: zhlou
 */

#ifndef EXPPARALLEL_H_
#define EXPPARALLEL_H_

#include "exponential.h"
#include "MPIState.h"

class expParallel
{
private:
    exponential serExp;
    const MPIState &mpi;
public:
    class Param
    {
    public:
        exponential::Param serParam;
        Param(xmlNode *root);
    };
    static const char *name;
    expParallel(const Param &param, const MPIState &mpiState);
    ~expParallel(){}
    void initStats(aState state){serExp.initStats(state);}
    void updateInitStep(bool accept, aState state) {serExp.updateInitStep(accept,state);}
    void resetSegmentStats(){serExp.resetSegmentStats();}
    void updateStep(bool accept, aState state){serExp.updateStep(accept,state);}
    double updateS(aState state){return serExp.updateS(state);}
    bool inSegment(aState state) {return serExp.inSegment(state);}
    void updateSegment(aState state){serExp.updateSegment(state);}
    void setDebug(debugStatus st, const char *outname=NULL){serExp.setDebug(st, outname);}
};

#endif /* EXPPARALLEL_H_ */
