#include "PWLE.h"
#include <cmath>

PWLE::PWLE(Param &param):
    nBins(param.nBins),
    eMin(param.eMin),
    binWidth(param.binWidth),
    syncFreq(syncFreq),
    uSteps(0)
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
    uSteps(orig.uSteps)
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

void PWLE::update(double eVal, double weight)
{
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
}
