#ifndef PWLE_H
#define PWLE_H

#include "MPIState.h"

class PWLE {
public:
    struct Param {
        unsigned int nBins;
        double eMin;
        double binWidth;
        MPIState &mpi;
        unsigned int syncFreq;
    };
    PWLE(Param &param);
    PWLE(const PWLE &orig);
    ~PWLE();

    void update(double eVal, double weight);
    double getV(double eVal) const;
private:
    const unsigned int nBins;
    const double eMin;
    const double binWidth;
    const unsigned int syncFreq;
    unsigned int uSteps;
    double *histLocal;
    double *histGlobal;
    int getBinNum(double eVal) const {return (int)((eVal-eMin)/binWidth);} 
}

#endif
