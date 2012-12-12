/*
 * invLinearFit.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: zhlou
 */

#include "invLinearFit.h"

invLinearFit::invLinearFit(double w) :
        weight(w)
{
    reset();


}

invLinearFit::~invLinearFit()
{
    // TODO Auto-generated destructor stub
}

void invLinearFit::reset()
{
    usx = 0.;
    usy = 0.;
    usxx = 0.;
    usxy = 0.;
    usyy = 0.;
    usum = 0.;
    A = 0.;
    B = 0.;
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
