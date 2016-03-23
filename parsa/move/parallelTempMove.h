#ifndef PARALLELTEMPMOVE_H
#define PARALLELTEMPMOVE_H

#include <mpi.h>
#include "mix/mixState.h"
#include "MPIState.h"
#include "move/tempFeedback.h"

template<class Problem>
class parallelTempMove : public tempFeedback<Problem>
{
public:
    parallelTempMove(Problem &in_problem, unirandom &in_rnd, const ptree &root,
            const MPIState &mpiState);
    ~parallelTempMove();
    int getWinner();
    void processMix(const mixState &ms, const aState &state)
    { if (ms.doesMix()) this->energy = state.energy;}
    static const char *name;
    void readState(const ptree &root); // it overrides the base class. bite me
protected:
    void collectMoveStats();
private:
    const MPIState &mpi;
};

#include "parallelTempMove.hpp"


#endif
