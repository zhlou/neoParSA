/* 
 * File:   DoS.h
 * Author: zhlou
 *
 * Created on July 10, 2015, 5:51 PM
 */

#ifndef DOS_H
#define	DOS_H

#include "unirandom.h"

template <class Problem, template<class> class Move, class Estimator>
class DoS
{
public:
    struct Param {
        long nSteps;
        double initWeight;
        typename Estimator::Param estParam;
    };
    DoS(Problem &problem, Move<Problem> &move, unirandom &in_rand, Param &param);
    void estimate();
    Estimator &getEstimator() {return estm;}
    void setWeight(double newWeight) {weight = newWeight;}
    double getWeight() const {return weight;}
private:
    Problem &problem;
    Move<Problem> &move;
    unirandom& rnd;
    Estimator estm;
    double weight;
    long nSteps;
    //double energy;
};

#include "DoS.hpp"

#endif	/* DOS_H */

