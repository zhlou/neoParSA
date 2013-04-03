/*
 * simpleAnnealer.h
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#ifndef SIMPLEANNEALER_H_
#define SIMPLEANNEALER_H_

#include "annealer.h"
#include "aState.h"
#include "dynDebug.h"
#include <libxml/tree.h>

class simpleSchedule
{
public:
    simpleSchedule(xmlNode *root);
    virtual ~simpleSchedule();
    void initStats(aState state);
    //double getInitS();
    //int getInitLoop();
    void updateInitStep(bool accept, aState state);
    virtual bool frozen(aState state);
    void resetSegmentStats(){};
    void updateStep(bool accept, aState state);
    double updateS(aState state);
    bool inSegment(aState state) {return !(state.step_cnt % max_rej);}
    void updateSegment(aState state);
    bool needMix(){return false;}
    void setDebug(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st,outname);}
    static const char* name;
protected:
    unsigned reject_cnt;
    double alpha;
    dynDebug debugOut;
    unsigned long max_rej;
    //int init_loop;
    //double init_S;
};

#endif /* SIMPLEANNEALER_H_ */
