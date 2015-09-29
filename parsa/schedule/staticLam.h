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
    std::vector<double> betaVec;
    std::vector<double> variance;
    dynDebug debugOut;
    const unsigned segLength;
    const double lambda;
    unsigned count;
    unsigned success;
    std::vector<double>::size_type size, i;
    double b0, bEnd, v0, cEnd;
    double alpha;

    double getVar(double beta);

public:
    class Param {
    public:
        debugStatus st;
        const char *logname;
        unsigned segLength;
        double lambda;
        char *filename;
        Param(xmlNode *root, debugStatus in_st=ignore, const char *logname=NULL);
    }
    staticLam(Param &param);
    ~staticLam();
    void initStats(const aState &){}
    void initStats(double, double, double, const aState&){}
    void updateInitStep(bool accept, const aState &state){updateStep(accept, state);}
    void resetSegmentStats(){}
    void updateStep(bool accept, const aState &state);
    double updateS(const aState &state);
    void updateSegment(const aState &state);
    void updateStats(const aState &state);
    void setDebug(debugStatus st, const char *outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
    static const char * name;
}



#endif
