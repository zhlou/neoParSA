/*
 *       Filename:  staticLam.h
 *        Created:  09/24/2015 03:47:26 PM
 *         Author:  Zhihao Lou
 *
 *                  This cooling schedule takes a text file with two columns
 *                  each line representing inverse temperature and
 *
 *           Note:  This cooling schedule relies on the properties that the
 *                  variance is constant at high temperature (small s or beta)
 *                  and is proportional to constant/beta^2 for low temperature
 *                  (large s).  This behavior is consistant with Gamma energy
 *                  density model and also theoretical values from the
 *                  Rastrigin functions.
 */

#ifndef STATICLAM_H
#define STATICLAM_H

#include <vector>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "dynDebug.h"
#include "aState.h"

class staticLam {
private:
    std::vector<double> betaVec;
    std::vector<double> variance;
    dynDebug debugOut;
    const unsigned segLength;
    const double lambda;
    const double minRate;
    unsigned count;
    unsigned success;
    std::vector<double>::size_type size, i;
    double b0, bEnd, v0, cEnd;
    double alpha;
    const int adjustAlpha;

    double getVar(double beta);
    void calcStats(unsigned nsteps, const aState &state);

public:
    class Param {
    public:
        debugStatus st;
        const char *logname;
        unsigned segLength;
        double lambda;
        double minRate;
        int adjustAlpha;
        const char *filename;
        Param(const ptree &root, debugStatus in_st=ignore, const char *logname=NULL);
    };
    staticLam(Param &param);
    ~staticLam(){}
    void initStats(const aState &);
    void initStats(double, double, double, const aState&);
    void updateInitStep(bool accept, const aState &state){updateStep(accept, state);}
    void updateStep(bool accept, const aState &state);
    double updateS(const aState &state);
    //void updateSegment(const aState &state);
    //void resetSegmentStats(){}
    void updateStats(const aState &state);
    void setDebug(debugStatus st, const char *outname=NULL)
    {
        debugOut.setDebug(st, outname);
    }
    static const char * name;
};



#endif
