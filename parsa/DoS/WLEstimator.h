/* 
 * File:   WLEsitmator.h
 * Author: zhlou
 *
 * Created on July 10, 2015, 2:15 PM
 */

#ifndef WLESTIMATOR_H
#define	WLESTIMATOR_H

#include <ostream>
#include <iostream>

class WLEstimator {
public:
    struct Param{
        unsigned int nBins;
        double eMin;
        double binWidth;
    };
    WLEstimator(Param &param);
    WLEstimator(const WLEstimator& orig);
    ~WLEstimator();
    void update(double eVal, double weight);

    double getV(double eVal) const;
    void printHist(std::ostream &out = std::cout) const;
    void saveHist(const char *filename) const;
    void readHist(const char *filename);
private:
    const unsigned int nBins;
    const double eMin;
    const double binWidth;
    double *hist;
    int getBinNum(double eVal) const {return (int)((eVal-eMin)/binWidth);} 

};

#endif	/* WLESITMATOR_H */

