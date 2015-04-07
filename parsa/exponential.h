/*
 * exponential.h
 *
 *  Created on: Dec 9, 2012
 *      Author: zhlou
 */

#ifndef EXPONENTIAL_H_
#define EXPONENTIAL_H_

#include "annealer.h"
#include "aState.h"
#include "dynDebug.h"
#include <libxml/tree.h>

class exponential
{
public:
    class Param{
    public:
        unsigned long max_rej;
        unsigned segLength;
        double alpha;
        debugStatus st;
        const char *outname;
        Param(xmlNode *root, debugStatus st=ignore, const char *outname=NULL);
    };
    exponential(const Param &param);
    // exponential(xmlNode *root);
    ~exponential();
    void initStats(aState state);
    //double getInitS();
    //int getInitLoop();
    void updateInitStep(bool accept, aState state);
    // virtual bool frozen(aState state);
    void resetSegmentStats(){};
    void updateStep(bool accept, aState state);
    double updateS(aState state);
    bool inSegment(aState state) {return (state.step_cnt % segLength);}
    void updateSegment(aState state);
    void updateStats(aState state);
    // bool needMix(){return false;}
    void setDebug(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st,outname);}
    static const char* name;
private:
    unsigned reject_cnt;
    double alpha;
    dynDebug debugOut;
    unsigned long max_rej;
    unsigned segLength;
    unsigned step_cnt;
    //int init_loop;
    //double init_S;
    // xmlNode *xmlsection;
};

#endif /* EXPONENTIAL_H_ */
