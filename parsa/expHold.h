/*
 * expHold.h
 *
 *  Created on: May 5, 2013
 *      Author: zhlou
 */

#ifndef EXPHOLD_H_
#define EXPHOLD_H_

#include "aState.h"
#include "dynDebug.h"
#include <libxml/tree.h>

class expHold {
private:
    double target_s;
    double alpha;
    dynDebug debugOut;
public:
    class Param {
    public:
        double target_s;
        double alpha;
        debugStatus st;
        const char *outname;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *name=NULL);
    };
    expHold(Param param) : target_s(param.target_s), alpha(param.alpha),
            debugOut(param.st, param.outname)
    {};
    ~expHold();
    void initStats(const aState &){}
    void updateInitStep(bool, const aState &){}
    void resetSegmentStats(){}
    void updateStep(bool, const aState &) {}
    double updateS(const aState &state);
    bool inSegment(aState) {return false;}
    void updateSegment(const aState &) {}
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
    static const char * name;
};

#endif /* EXPHOLD_H_ */
