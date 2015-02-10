/* 
 * File:   moveControlCore.h
 * Author: zhlou
 *
 * Created on February 5, 2015, 2:06 PM
 */

#ifndef MOVECONTROLCORE_H
#define	MOVECONTROLCORE_H

#include "MPIState.h"
#include "unirandom.h"
#include "dynDebug.h"

class moveControlCore
{
private:
    const int nparams;
    const double gain;
    const double target;
    
    long *success;
    long *moves;
    double *actualThetas;
    double *thetaBars;
    double *thetaMins;
    double *thetaMaxs;
    double varTheta;
    
    const MPIState &mpi;
    unirandom &rnd;
    dynDebug &debugOut;
    void setActualTheta(int index);
public:
    moveControlCore(int nparams, double gain, double target, double initTheta,
            double thetaMin, double thetaMax, const MPIState &mpi, 
            unirandom &rnd, dynDebug &debugOut);
    ~moveControlCore();
    double genMove(int index) {
        ++moves[index];
        return rnd.laplace(actualThetas[index]);
    }
    void accept(int index) {++ success[index];};
    void reject(int index) {};
    void setVarTheta(double newVarTheta);
    double getVarTheta() const {return varTheta;}
    void setThetaBar(int index, double thetaBar);
    double getThetaBar(int index) const {return thetaBars[index];}
    double getActualTheta(int index) const {return actualThetas[index];}
    void moveControl();
    
};

#endif	/* MOVECONTROLCORE_H */

