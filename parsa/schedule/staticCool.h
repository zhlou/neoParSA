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
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "dynDebug.h"
#include "aState.h"




class staticCool {
private:
    std::vector<double> schedule;
    dynDebug debugOut;
    const unsigned segLength;
    unsigned step_cnt;
    std::vector<double>::size_type size, i;
public:
    class Param {
    public:
        const char *scheduleName;
        const char *outname;
        debugStatus st;
        unsigned segLength;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *logname=NULL);
        Param(const ptree &root, debugStatus in_st=ignore, const char *logname=NULL);
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
