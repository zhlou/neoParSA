/* 
 * File:   WLEsitmator.cpp
 * Author: zhlou
 * 
 * Created on July 10, 2015, 2:15 PM
 */

#include "WLEstimator.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>

WLEstimator::WLEstimator(Param &param) :
        nBins(param.nBins),
        eMin(param.eMin),
        binWidth(param.binWidth)
{
    hist = new double[nBins];
    for (int i = 0; i < nBins; ++i)
        hist[i] = 0;
}

WLEstimator::WLEstimator(const WLEstimator& orig) :
        nBins(orig.nBins),
        eMin(orig.eMin),
        binWidth(orig.binWidth),
        hist(orig.hist)
{ 
    hist = new double[nBins];
    for (int i = 0; i < nBins; ++i)
        hist[i] = orig.hist[i];
}

WLEstimator::~WLEstimator() 
{ 
    delete []hist;
}

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


void WLEstimator::saveHist(const char *filename) const
{
    std::ofstream fp(filename, std::ios::binary);
    fp.write((const char*)&hist[0], nBins*sizeof(double));
    fp.close();
}

void WLEstimator::readHist(const char *filename)
{
    std::ifstream fp(filename, std::ios::binary);
    fp.read((char*)&hist[0], nBins*sizeof(double));
    if (!fp) {
        throw std::runtime_error("reading hist from file failed");
    }
    fp.close();
}
