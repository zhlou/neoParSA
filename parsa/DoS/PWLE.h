#ifndef PWLE_H
#define PWLE_H

#include "MPIState.h"
#include <iostream>

class PWLE {
public:
    struct Param {
        unsigned int nBins;
        double eMin;
        double binWidth;
        MPIState *mpi;
        unsigned int syncFreq;
    };
    PWLE(Param &param);
    PWLE(const PWLE &orig);
    ~PWLE();

    void update(double eVal, double weight);
    double getV(double eVal) const;
    void printHist(std::ostream &out = std::cout) const;
    void saveHist(const char *filename) const;
    void readHist(const char *filename);
    void syncHist();
private:
    const unsigned int nBins;
    const double eMin;
    const double binWidth;
    const unsigned int syncFreq;
    MPIState &mpi;
    unsigned int uSteps;
    double *histLocal;
    double *histGlobal;
    bool localModified;
    int getBinNum(double eVal) const {return (int)((eVal-eMin)/binWidth);} 

};

#endif
