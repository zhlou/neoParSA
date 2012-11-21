#ifndef TSP_PROBLEM_H
#define TSP_PROBLEM_H

#include "movable.h"
#include "tsp.h"
class tsp_reorder: public abstract_param {
    public:
        void generate_tweak(double theta_bar);
        void restore_tweak();
        tsp_reorder(tsp *from_problem);
        ~tsp_reorder();
    private:
        tsp *thetsp;
        unsigned seed;
        int ncities;
};

class tsp_problem: public movable {
    public:
        double get_score();
        tsp_problem(tsp *ext_tsp);
        ~tsp_problem();
    private:
        tsp *thetsp;

};

#endif

