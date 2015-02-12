/* 
 * File:   FBMoveFixedMix.h
 * Author: zhlou
 *
 * Created on January 29, 2015, 3:48 PM
 */

#ifndef FBMOVEINTERVALMIX_H
#define	FBMOVEINTERVALMIX_H

#include <libxml/tree.h>
#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include "Mixing.h"
#include "unirandom.h"
#include "moveControlCore.h"

template<class Problem>
class FBMoveIntervalMix {
public:

    class Param {
    public:
        int mix_interval;
        int move_interval;
        double move_gain;
        double target;
        double initTheta;
        double thetaMin;
        double thetaMax;
        double mix_target;
        double varConst;
        Param();
        Param(xmlNode *root);

    };
    FBMoveIntervalMix(Problem &in_problem, const MPIState &mpiState,
            unirandom &in_rnd, const Param &param);
    ~FBMoveIntervalMix();
    double get_score() { return energy;}
    double propose();
    void accept();
    void reject();
    mixState Mix(aState &state);
    void setDebug(debugStatus st, const char* outname=NULL) {
        debugOut.setDebug(st, outname);}
    void setMixLog(debugStatus st, const char* outname=NULL) {
        mixLog.setDebug(st, outname);
    }
    static const char * name;
private:
    Problem &problem;
    unirandom& rnd;
    const MPIState &mpi;
    dynDebug debugOut;
    dynDebug mixLog;
    const int nparams;
    const int mix_interval;
    const int move_interval;

    int index;
    long sweep;
    double energy;
    double prev_energy;

    int tau_count;
    
    Mixing<Problem> mix;
    moveControlCore moveCore;
    double mix_target;
    double varConst;
    int *adoptArray;

    void move_control();

};


#include "FBMoveIntervalMix.hpp"

#endif	/* FBMOVEFIXEDMIX_H */

