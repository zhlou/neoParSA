// created 2/16/2016
#ifndef SCHEDULETHETAMOVEP_H
#define SCHEDULETHETAMOVEP_H

#include "scheduleThetaMove.h"
#include "MPIState.h"

template <class Problem>
class scheduleThetaMoveP : public scheduleThetaMove<Problem>
{
public:
    scheduleThetaMoveP(Problem &in_problem, unirandom &in_rnd, const ptree &root, MPIState &mpi) : scheduleThetaMove<Problem>(in_problem, in_rnd, root) {}
};


#endif
