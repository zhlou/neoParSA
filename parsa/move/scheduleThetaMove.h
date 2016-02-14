// Created 2/13/2016
#ifndef SCHEDULETHETAMOVE_H
#define SCHEDULETHETAMOVE_H

#include <vector>

#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "dynDebug.h"
#include "unirandom.h"
#include "aState.h"

template <class Problem>
class scheduleThetaMove
{
public:
    scheduleThetaMove(Problem &in_problem, unirandom &in_rnd, const ptree &root);
    virtual ~scheduleThetaMove();
    double get_score(){return energy;}
    double propose(double s);
    void accept();
    void reject();
    void setDebug(debugStatus st, const char*outname=NULL)
    {debugOut.setDebug(st, outname);}
    void writeState(ptree &) const {}
    void readState(const ptree &root) {}
    double forceUpdateEnergy() {return energy = problem.get_score();}
    static const char *name;
protected:
    int nparams;
    Problem &problem;
    dynDebug debugOut;
    double energy;
    std::vector<unsigned> moves;
    std::vector<unsigned> success;
    //thetaTab[][0] is inverse temp and thetaTab[][1:nparams+1] are theta_bars;
    std::vector<std::vector<double> > thetaTab;
private:
    unsigned index;
    unsigned currentTabRow;
    unirandom &rnd;
    double prev_energy;
    unsigned interval;
    unsigned long sweep;
    double getThetaByS(double s, unsigned index);
};

#include "scheduleThetaMove.hpp"

#endif
