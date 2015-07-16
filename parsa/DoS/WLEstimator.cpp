/* 
 * File:   WLEsitmator.cpp
 * Author: zhlou
 * 
 * Created on July 10, 2015, 2:15 PM
 */

#include "WLEstimator.h"
#include <cmath>
#include <iostream>

WLEstimator::WLEstimator(Param &param) :
        nBins(param.nBins),
        eMin(param.eMin),
        binWidth(param.binWidth),
        hist(nBins,0)
{
}

WLEstimator::WLEstimator(const WLEstimator& orig) :
        nBins(orig.nBins),
        eMin(orig.eMin),
        binWidth(orig.binWidth),
        hist(orig.hist)
{ 
}

WLEstimator::~WLEstimator() { }

void WLEstimator::update(double eVal, double weight)
{
    int binN = getBinNum(eVal);
    if ( (unsigned) binN < nBins )
        hist[binN] += weight;
}

double WLEstimator::getV(double eVal) const
{
    int binN = getBinNum(eVal);
    if ((unsigned) binN < nBins) {
        return hist[binN];
    } else {
        return HUGE_VAL;
    }
}

void WLEstimator::printHist(std::ostream &out) const 
{
    for (int i = 0; i < nBins; ++i) {
        out << eMin+i*binWidth << "    " << hist[i] << std::endl;
    }
}
