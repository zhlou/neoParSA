/*
 * simPLam.cpp
 *
 *  Created on: Dec 17, 2013
 *      Author: zhlou
 */

#include <cmath>
#include "simPLam.h"

const char * simPLam::name = "simPLam";

simPLam::Param::Param(xmlNode* root, debugStatus in_st, const char* name) {
}

void simPLam::updateStep(bool accept, const aState &state)
{
    if (accept)
        ++ success;
    sum += state.energy;
    sumsq += (state.energy * state.energy);
}

simPLam::simPLam(Param param, const MPIState&) {
}

simPLam::~simPLam() {
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
}

double simPLam::updateS(const aState& state)
{
    return state.s + lambda * alpha / (sd * sd * sd * state.s);
}

void simPLam::calcStats(int nsteps)
{
    double from[2] = { sum, sumsq }, to[2];
    double N = (double) (mpiState.nnodes) * nsteps;
    MPI_Allreduce(from, to, 2, MPI_DOUBLE, MPI_SUM, mpiState.comm);
    mean = to[0] / N;
    sd = std::sqrt((to[1] - N * to[0] * to[0]) / N);
    acc_ratio = (double) (success) / N;
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}

void simPLam::updateSegment(const aState&)
{
    calcStats(proc_tau);
}

void simPLam::setDebug(debugStatus st, const char* outname) {
}
