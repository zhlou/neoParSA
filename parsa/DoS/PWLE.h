#ifndef PWLE_H
#define PWLE_H

#include "MPIState.h"
#include "WLEstimator.h"

class PWLE : public WLEstimator {
public:
    struct Param {
        WLEstimator::Param serialParam;
        MPIState &mpi;
        unsigned int syncFreq;
    };
    PWLE(Param &param);
    PWLE(const PWLE &orig);
    ~PWLE();
}

#endif
