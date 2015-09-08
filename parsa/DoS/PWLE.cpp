#include "PWLE.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>

PWLE::PWLE(Param &param):
    nBins(param.nBins),
    eMin(param.eMin),
    binWidth(param.binWidth),
    syncFreq(syncFreq),
    mpi(*param.mpi),
    uSteps(0),
    localModified(false)
{
    histLocal = new double[nBins];
    histGlobal = new double[nBins];
    for (int i = 0; i < nBins; ++i) {
        histLocal[i] = 0;
        histGlobal[i] = 0;
    }
}

PWLE::PWLE(const PWLE &orig) : 
    nBins(orig.nBins),
    eMin(orig.eMin),
    binWidth(orig.binWidth),
    syncFreq(orig.syncFreq),
    mpi(orig.mpi),
    uSteps(orig.uSteps),
    localModified(orig.localModified)
{
    histLocal = new double[nBins];
    histGlobal = new double[nBins];
    for (int i = 0; i < nBins; ++i) {
        histLocal[i] = orig.histLocal[i];
        histGlobal[i] = orig.histGlobal[i];
    }
}

PWLE::~PWLE()
{
    delete []histGlobal;
    delete []histLocal;
}

void PWLE::syncHist()
{
    MPI_Allreduce(MPI_IN_PLACE, histLocal, nBins, MPI_DOUBLE, MPI_SUM, mpi.comm);
    for (int i = 0; i < nBins; ++i) {
        histGlobal[i] += histLocal[i];
        histLocal[i] = 0;
    }
    localModified = false;
}

void PWLE::update(double eVal, double weight)
{
    int binN = getBinNum(eVal);
    if ( (unsigned) binN < nBins )
        histLocal[binN] += weight;
    localModified = true;
    ++ uSteps;
    if (syncFreq == uSteps) {
        uSteps = 0;
        syncHist();
    }
}

double PWLE::getV(double eVal) const
{
    int binN = getBinNum(eVal);
    if ((unsigned) binN < nBins) {
        return (histLocal[binN] + histGlobal[binN]);
    } else {
        return HUGE_VAL;
    }
}

void PWLE::readHist(const char* filename)
{
    std::ifstream fp(filename, std::ios::binary);
    fp.read((char*)&histGlobal[0], nBins*sizeof(double));
    if (!fp) {
        throw std::runtime_error("reading hist from file failed");
    }
    fp.close();
    for (int i = 0; i < nBins; ++i) {
        histLocal[i] = 0;
    }
    localModified = false;
}

void PWLE::saveHist(const char* filename) const
{
    if (localModified)
        std::cerr << "Warning: data in histLocal will be ignored" << std::endl;
    std::ofstream fp(filename, std::ios::binary);
    fp.write((const char*)&histGlobal[0], nBins*sizeof(double));
    fp.close();
}

void PWLE::printHist(std::ostream& out) const
{
    if (localModified)
        std::cerr << "Warning: data in histLocal will be ignored" << std::endl;
    for (int i = 0; i < nBins; ++i) {
        out << eMin+i*binWidth << "    " << histGlobal[i] << std::endl;
    }
}