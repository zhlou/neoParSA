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
#include "onePassMeanVar.h"
#include <libxml/tree.h>

class expHold {
private:
    unsigned segLength;
    double target_s;
    double alpha;
    dynDebug debugOut;
    onePassMeanVar energyStat;
    unsigned step_cnt;
public:
    class Param {
    public:
        unsigned segLength;
        double target_s;
        double alpha;
        debugStatus st;
        const char *outname;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *name=NULL);
    };
    expHold(Param param) : segLength(param.segLength), target_s(param.target_s),
            alpha(param.alpha), debugOut(param.st, param.outname), energyStat(),
            step_cnt(0)
    {};
    virtual ~expHold();
    void initStats(const aState &){}
    void initStats(double, double, double, const aState &){}
    void updateInitStep(bool, const aState &){}
    void resetSegmentStats(){energyStat.reset();}
    void updateStep(bool, const aState &state) {energyStat.update(state.energy);}
    double updateS(const aState &state);
    // bool inSegment(aState state) {return (state.step_cnt % segLength);}
    void updateSegment(const aState &state);
    void updateStats(const aState &state);
    void setDebug(debugStatus st, const char* outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
    static const char * name;
};

#endif /* EXPHOLD_H_ */
