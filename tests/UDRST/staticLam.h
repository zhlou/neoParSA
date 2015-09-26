/*
 *       Filename:  staticLam.h
 *        Created:  09/24/2015 03:47:26 PM
 *         Author:  Zhihao Lou
 *
 *           Note:  This cooling schedule relies heavily on the properties
 *                  of Rastrigin function so it has to be here.
 */

#ifndef STATICLAM_H
#define STATICLAM_H

#include <vector>
#include <libxml/tree.h>

#include "dynDebug.h"
#include "aState.h"

class staticLam {
private:
    std::vector<double> energy;
    std::vector<double> variance;
    dynDebug debugOut;
    const unsigned segLength;
    unsigned step_cnt;
    std::vector<double>::size_type size, i;
    const double lambda;
    double getVar(double beta);

public:
    class Param {
    public:
        char *filename;
        const char *logname;
        debugStatus st;
        unsigned segLength;
        double lambda;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *logname=NULL);
    }
    staticLam(Param &param);
    ~staticLam();
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
}



#endif
