/*
 * invLinearFit.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: zhlou
 */

#include "invLinearFit.h"

invLinearFit::invLinearFit(double w, double d0, double s0, double inA) :
        weight(w), A(inA)
{
    B = (1.0 / d0) - (A * s0);
    usum = 1.0;
    usx = 0.;
    usy = 1.0 / d0;
    usxx = 0.;
    usxy = 0.;
    usyy = usy * usy;
}

invLinearFit::~invLinearFit()
{
    // TODO Auto-generated destructor stub
}

void invLinearFit::reset()
{
    usx = 0.;
    usy = 1.;
    usxx = 0.;
    usxy = 0.;
    usyy = 1.;
    usum = 1.;
    A = 1.; //these two initialize to 1 so it always has value
    B = 1.;
}

void invLinearFit::fullUpdate(double d, double s)
{
    decay();
    partialUpdate(1.0, d, s);
    finishUpdate();

}

void invLinearFit::decay()
{
    usx *= weight;
    usy *= weight;
    usxx *= weight;
    usxy *= weight;
    usyy *= weight;
    usum *= weight;
}

void invLinearFit::partialUpdate(double p_weight, double d, double s)
{
    usyy += p_weight*d*d;
    usxy += p_weight*s*d;
    usy += p_weight*d;
    usx += p_weight*s;
    usxx += p_weight*s*s;
    usum += p_weight;
}

void invLinearFit::finishUpdate()
{
    A = (usum * usxy - usx * usy) / (usum * usxx - usx * usx);
    B = (usy - A * usx) / usum;
}

double invLinearFit::getEstimate(double s)
{
    return 1.0 / (A * s + B);
}
