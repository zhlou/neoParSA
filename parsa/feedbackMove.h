#ifndef MOVABLE_H
#define MOVABLE_H
#include <libxml/parser.h>

/*
 * This is the feedback move control template. It features perturbation of
 * single parameter at each move and proportion feedback control at the end
 * of each sweep to adjust the theta_bars so the acceptance ratio to the
 * reference value of 0.44. It needs to take the Problem class to have the
 * following interfaces
 *
 * class Problem
 * {
 * public:
 *      int getDimension(); // returns the dimension of the problem
 *      double get_score(); // returns the energy (score) of the problem
 *      void generateMove(int index, double theta_bar);
 *      // make the perturbation on parameter index with theta_bar
 *      void restoreMove(int index);
 * };
 *
 */
template <class Problem, class Debug>
class feedbackMove : public Debug
{
public:
    feedbackMove(Problem &in_problem, xmlNode *root=NULL);
    virtual ~feedbackMove();
    double get_score();
    double propose();
    void accept();
    void reject();
    virtual void doMix(){}; // do nothing in base class

protected:
    int nparams;
    Problem &problem;
    virtual void collectMoveStats(){}; // do nothing in base class
    void move_control();
    long *success;
    long *moves;
private:
    int index;

    double *theta_bars;
    double prev_energy;
    double energy;
    double move_gain;
    long sweep;
    int move_interval;
    static const double theta_min;
};



#include "feedbackMove.hpp"

#endif
