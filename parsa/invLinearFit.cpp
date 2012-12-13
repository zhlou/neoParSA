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

void invLinearFit::update(double d, double s)
{
    usx *= weight;
    usy *= weight;
    usxx *= weight;
    usxy *= weight;
    usyy *= weight;
    usum *= weight;

    usyy += d*d;
    usxy += s*d;
    usy += d;
    usx += s;
    usxx += s*s;
    usum += 1.0;

    A = (usum * usxy - usx * usy) / (usum * usxx - usx * usx);
    B = (usy - A * usx) / usum;
}

double invLinearFit::getEstimate(double s)
{
    return 1.0 / (A * s + B);
}
