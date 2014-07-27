/*
 * fixedTheta.h
 *
 *  Created on: Jul 26, 2014
 *      Author: zhlou
 */

#ifndef FIXEDTHETAMOVE_H_
#define FIXEDTHETAMOVE_H_

#include <libxml/parser.h>

#include "unirandom.h"
#include "dynDebug.h"

template <class Problem>
class fixedThetaMove
{
public:
    fixedThetaMove(Problem &in_problem, unirandom& in_rnd, xmlNode *root);
    ~fixedThetaMove();
    static const char *name;
    double get_score(){return energy;}
    double propose();
    void accept();
    void reject();
    void setDebug(debugStatus st, const char *outname=NULL)
    {debugOut.setDebug(st, outname);}
private:
    const int nparams;
    int index;
    unirandom& rnd;
    Problem &problem;
    dynDebug debugOut;
    double target;
    int logInterval;
    int counter;
    double energy;
    double prev_energy;
    double accumTheta;
    int accumAccept;
    void procStats();


};

#include "fixedThetaMove.hpp"

#endif /* FIXEDTHETAMOVE_H_ */
