/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"
#include <cmath>
using namespace std;

lam::lam(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    proc_tau = 100; // TODO for now. should be an input from xml later
    fit_mean = new invLinearFit(0.995);
    fit_sd = new invLinearFit(0.99);
    resetSegmentStats();
}

lam::~lam()
{
    delete fit_mean;
    delete fit_sd;
}
/*
 double lam::loop()
 {
 long unsigned step_cnt = 0;
 int accept = 0, i;
 double delta, vari;
 while (!frozen()) {
 accept = 0;
 vari = 0.;
 for (i = 0; i < proc_tau; i++) {
 if ((delta = move()) != 0.)
 accept ++;

 }

 }

 return problem->get_score();
 }
 */

bool lam::frozen()
{
    return false; //TODO need to add frozen condition
}

void lam::updateStep(bool is_accept, double delta)
{
    if (is_accept)
        success++;
    vari += energy - estimate_mean;
    mean += energy;
}

void lam::updateS()
{
}

void lam::resetSegmentStats()
{
    acc_ratio = 0.;
    success = 0;
    vari = 0.;
    mean = 0.;
}

bool lam::inSegment()
{
    if ((step_cnt % proc_tau) == 0) {
        return false;
    } else {
        return true;
    }
}

void lam::updateSegment()
{
    mean /= proc_tau;
    vari /= proc_tau;
    acc_ratio = (double)success / proc_tau;
    updateLam();
}

void lam::resetLam()
{
    fit_mean->reset();
    fit_sd->reset();
}

void lam::updateLam()
{
    fit_mean->update(1.0/mean, s);
    fit_sd->update(1.0/sqrt(vari),s);
    double d = (1.0 - acc_ratio) / (2.0 - acc_ratio);
    alpha = 4.0 * acc_ratio * d * d;
}
