/*
 * simPLam.cpp
 *
 *  Created on: Dec 17, 2013
 *      Author: zhlou
 */

#include <cmath>
#include <exception>
#include <stdexcept>
#include "simPLam.h"
#include "utils.h"


const char * simPLam::name = "simPLam";

simPLam::Param::Param(xmlNode* root, debugStatus in_st, const char* name) :
        st(in_st),
        outname(name)
{
    xmlNode *xmlsection = getSectionByName(root, "simPLam");
    if (xmlsection == NULL)
        throw std::runtime_error(std::string("Error: fail to find section simPLam"));
    proc_tau = getPropInt(xmlsection,"tau");
    lambda = getPropDouble(xmlsection, "lambda");

}

void simPLam::updateStep(bool accept, const aState &state)
{
    if (accept)
        ++ success;
    sum += state.energy;
    sumsq += (state.energy * state.energy);
}

simPLam::simPLam(Param param, const MPIState& mpi) :
        proc_tau(param.proc_tau),
        lambda(param.lambda),
        mpiState(mpi),
        acc_ratio(0.),
        alpha(0.),
        mean(0.),
        sd(0.),
        sum(0.),
        sumsq(0.),
        success(0)
{

}

simPLam::~simPLam()
{
}

void simPLam::initStats(const aState& state)
{
    calcStats(state.step_cnt);
}

void simPLam::updateInitStep(bool accept, const aState &state)
{
    updateStep(accept, state);
}

void simPLam::resetSegmentStats()
{
    sum = 0;
    sumsq = 0;
    success = 0;
}

double simPLam::updateS(const aState& state)
{
    return state.s + lambda * alpha / (sd * sd * sd * state.s * state.s);
}

void simPLam::calcStats(int nsteps)
{
    double from[3] = { sum, sumsq, success }, to[3];
    double N = (double) (mpiState.nnodes) * nsteps;
    MPI_Allreduce(from, to, 3, MPI_DOUBLE, MPI_SUM, mpiState.comm);
    mean = to[0] / N;
    sd = std::sqrt((to[1] - mean * to[0]) / N);
    acc_ratio = (to[2]) / N;
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}

void simPLam::updateSegment(const aState &state)
{
    calcStats(proc_tau);
    debugOut << state.step_cnt << " " << state.s << " "
             << lambda * alpha / (state.s * state.s * sd * sd * sd)
             << " " << mean << " " << sd << " " << acc_ratio << " "
             << alpha << std::endl;
}
