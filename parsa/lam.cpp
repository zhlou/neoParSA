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
    acc_ratio = 0.;
    accept = 0;
    nvari = 0.;

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
}

void lam::updateStats(bool is_accept, double delta)
{
    if (is_accept)
        accept ++;
    nvari += delta * delta;
    if (step_cnt % proc_tau == 0) {

    }

}

void lam::cool_s()
{
}
