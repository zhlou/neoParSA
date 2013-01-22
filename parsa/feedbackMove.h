#ifndef MOVABLE_H
#define MOVABLE_H
#include <libxml/parser.h>
/*
class abstract_param {
    public: 
        virtual void generate_tweak(double theta_bar) = 0;
        virtual void restore_tweak() = 0;
        virtual ~abstract_param(); // need this on obj w/ virtual func.
};
*/

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
template <class Problem>
class feedbackMove{
public:
    feedbackMove(Problem &in_problem, xmlNode *root=NULL);
    virtual ~feedbackMove();
    double get_score();
    double propose();
    void accept();
    void reject();

protected:
    int nparams;
    Problem &problem;
    virtual void move_control();
private:
    int index;
    long *success;
    long *moves;
    double *theta_bars;
    double prev_energy;
    double energy;
    double move_gain;
    long sweep;
    int move_interval;
    static const double theta_min;
};

/*
class movable {
    public:
        double propose_move();
        void accept_move();
        void reject_move();
        void init_stats();
        void set_theta(int id, double theta);
        virtual double get_score() = 0;
        virtual ~movable();
        movable(int np, xmlNode *root=NULL);

    protected:
        int nparams;
        //Problem &problem;
        abstract_param **params;

    private:
        void init(int np);
        void move_control();
        int index;
        long *success;
        long *moves;
        double *theta_bars;
        double energy;
        double prev_eng;
        double move_gain;
        long sweep;
        int move_interval;
        static const double theta_min;
};
*/


#include "feedbackMove.hpp"

#endif
