/* 
 * File:   staticCool.h
 * Author: zhlou
 *
 * Created on September 21, 2015, 10:00 AM
 */

#ifndef STATICCOOL_H
#define	STATICCOOL_H

#include <vector>

#include <libxml/tree.h>

#include "dynDebug.h"
#include "aState.h"
#include "utils/vectorUtils.h"




class staticCool {
private:
    std::vector<double> schedule;
    dynDebug debugOut;
public:
    class Param {
    public:
        char *scheduleName;
        const char *outname;
        debugStatus st;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *logname=NULL);
    };
    staticCool(Param &param);
    ~staticCool();
    void initStats(const aState &){}
    void initStats(double, double, double, const aState&){}
    void updateInitStep(bool, const aState &){}
    void resetSegmentStats(){}
    void updateStep(bool, const aState &state){}
    double updateS(const aState &state);
    void updateSegment(const aState &state);
    void updateStats(const aState &state);
    void  setDebug(debugStatus st, const char *outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
    static const char * name;
    
    
    
};

#endif	/* STATICCOOL_H */

