#include "moveControlCore.h"

#include <cmath>

moveControlCore::moveControlCore(int nparams, double gain, double target, 
        double initTheta, double thetaMin, double thetaMax, 
        const MPIState& mpi, unirandom &rnd, dynDebug &debugOut) : nparams(nparams), gain(gain), 
        target(target), varTheta(0), mpi(mpi), rnd(rnd), debugOut(debugOut)
{
    success = new long[nparams];
    moves = new long[nparams];
    actualThetas = new double[nparams];
    thetaBars = new double[nparams];
    thetaMins = new double[nparams];
    thetaMaxs = new double[nparams];
    for (int i = 0; i < nparams; ++i) {
        success[i] = 0;
        moves[i] = 0;
        actualThetas[i] = initTheta;
        thetaBars[i] = initTheta;
        thetaMins[i] = thetaMin;
        thetaMaxs[i] = thetaMax;
    }
    
}

moveControlCore::~moveControlCore() 
{
    delete[] thetaMaxs;
    delete[] thetaMins;
    delete[] thetaBars;
    delete[] actualThetas;
    delete[] moves;
    delete[] success;
}

void moveControlCore::setThetaBar(int index, double thetaBar)
{
    thetaBars[index] = thetaBar;
    setActualTheta(index);
}

void moveControlCore::setVarTheta(double newVarTheta) {

    varTheta = newVarTheta;
    for (int i = 0; i < nparams; ++ i) {
        setActualTheta(i);
    }
}


void moveControlCore::moveControl()
{
    MPI_Allreduce(MPI_IN_PLACE, success, nparams, MPI_LONG, MPI_SUM, mpi.comm);
    MPI_Allreduce(MPI_IN_PLACE, moves, nparams, MPI_LONG, MPI_SUM, mpi.comm);

    for (int i = 0; i < nparams; ++i) {
        double accRatio = (double) success[i] / (double) moves[i];
        double x = std::log(thetaBars[i]);
        debugOut << "\t" << accRatio;
        x += gain * (accRatio - target);
        setThetaBar(i,exp(x));
        success[i] = 0;
        moves[i] = 0;
    }
}

void moveControlCore::setActualTheta(int index)
{
    if (0. == varTheta)
        actualThetas[index] = thetaBars[index];
    else
        actualThetas[index] = std::exp(rnd.randn(std::log(thetaBars[index]), varTheta));
    if (actualThetas[index] < thetaMins[index])
        actualThetas[index] = thetaMins[index];
    if (actualThetas[index] > thetaMaxs[index])
        actualThetas[index] = thetaMaxs[index];
}
