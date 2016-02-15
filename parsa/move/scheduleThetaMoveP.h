// created 2/16/2016
#ifndef SCHEDULETHETAMOVEP_H
#define SCHEDULETHETAMOVEP_H

#include "scheduleThetaMove.h"
#include "MPIState.h"
#include "mix/mixState.h"

template <class Problem>
class scheduleThetaMoveP : public scheduleThetaMove<Problem>
{
public:
    scheduleThetaMoveP(Problem &in_problem, unirandom &in_rnd, const ptree &root, const MPIState &mpi) : scheduleThetaMove<Problem>(in_problem, in_rnd, root) {}
    void processMix(const mixState &ms, const aState &state)
    {if (ms.doesMix()) this->energy = state.energy;}
};


#endif
