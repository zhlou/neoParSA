/*
 * lam.cpp
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#include "lam.h"

lam::lam(movable *theproblem, xmlNode *root) :
        annealer(theproblem, root)
{
    proc_tau = 100; // TODO for now. should be an input from xml later
    resetStats();
}

lam::~lam()
{
    // TODO Auto-generated destructor stub
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
        accept++;
    vari += delta * delta;
    mean += delta;
}

void lam::cool_s()
{
}

void lam::resetStats()
{
    acc_ratio = 0.;
    accept = 0;
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
}
