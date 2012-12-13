/*
 * invLinearFit.h
 *
 *  Created on: Dec 11, 2012
 *      Author: zhlou
 */

#ifndef INVLINEARFIT_H_
#define INVLINEARFIT_H_

class invLinearFit
{
public:
    invLinearFit(double w, double d0, double s0, double inA);
    virtual ~invLinearFit();
    void reset();
    void update(double d, double s);
    double getEstimate(double s);
private:
    double weight;
    double usx, usy, usxx, usxy, usyy, usum;
    double A, B;
};

#endif /* INVLINEARFIT_H_ */
