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
    expHold(xmlNode *root);
    ~expHold();
    void initStates(const aState &){}
    void updateInitSte(bool, const aState &){}
    void resetSegmentStats(){}
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
