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
    w_a = 0.995;
    w_b = 0.99;
    resetSegmentStats();
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
}

void lam::resetLam()
{
    usx = 0.;
    usy = 0.;
    usxx = 0.;
    usxy = 0.;
    usyy = 0.;
    usum = 0.;

    vsx = 0.;
    vsy = 0.;
    vsxx = 0.;
    vsxy = 0.;
    vsyy = 0.;
    vsum = 0.;

}

void lam::updateLam()
{
    double d = 1.0/mean;
    usx *= w_a;
    usy *= w_a;
    usxx *= w_a;
    usxy *= w_a;
    usyy *= w_a;
    usum *= w_a;

    usyy += d*d;
    usxy += s*d;
    usy += d;
    usx += s;
    usxx += s*s;
    usum += 1.0;

    d = 1.0 / sqrt(vari);
    vsx *= w_b;
    vsy *= w_b;
    vsxx *= w_b;
    vsxy *= w_b;
    vsyy *= w_b;
    vsum *= w_b;

    vsyy += d*d;
    vsxy += s*d;
    vsy += d;
    vsx += s;
    vsxx += s*s;
    vsum += 1.0;


}
