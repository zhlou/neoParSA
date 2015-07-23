/* 
 * File:   DoS.hpp
 * Author: zhlou
 *
 * Created on July 10, 2015, 6:06 PM
 */

#ifndef DOS_HPP
#define	DOS_HPP

#include <cmath>

template <class Problem, template<class> class Move, class Estimator>
DoS<Problem, Move, Estimator>::DoS(Problem& problem, Move<Problem>& move, unirandom& in_rand, Param& param) : 
        problem(problem),
        move(move),
        rnd(in_rand),
        estm(param.estParam),
        weight(param.initWeight),
        nSteps(param.nSteps)
{
    //energy = move.get_score();
}

template <class Problem, template<class> class Move, class Estimator>
void DoS<Problem, Move, Estimator>::estimate()
{
    double energy = move.get_score();
    double Vold = estm.getV(energy);
    double Vnew = Vold;
    double energyNew = energy;
    double crit = 0.0, delta = 0.0;
    for (long i = 0; i < nSteps; ++i) {
        energyNew = energy + move.propose();
        Vold = estm.getV(energy);
        Vnew = estm.getV(energyNew);
        delta = Vold - Vnew;
        crit = std::exp(delta);
        if ( rnd.random() < crit) {
            move.accept();
            energy = energyNew;
            Vold = Vnew;
        } else {
            move.reject();
        }
        estm.update(energy, weight);
    }
}

#endif	/* DOS_HPP */

