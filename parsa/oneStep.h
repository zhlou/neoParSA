/*
 * oneStep.h
 *
 *  Created on: Apr 22, 2013
 *      Author: zhlou
 */

#ifndef ONESTEP_H_
#define ONESTEP_H_

#include "aState.h"
#include "dynDebug.h"
#include <libxml/tree.h>

class oneStep {
private:
    double target_s;
    dynDebug debugOut;
public:
    class Param {
    public:
        double target_s;
        debugStatus st;
        const char *outname;
        Param(xmlNode *root, debugStatus st=ignore, const char *outname=NULL);
    };
    oneStep(Param param);
    ~oneStep();
    void initStats(aState &state) {state.s = target_s;}
    void updateInitStep(bool, aState){}
    void resetSegmentStats(){}
    void updateStep(bool, aState){}
    double updateS(const aState &state){return state.s;}
    bool inSegment(aState) {return false;}
    void updateStats(aState){}
    void setDebug(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st,outname);}
    static const char * name;


};

#endif /* ONESTEP_H_ */
