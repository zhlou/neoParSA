#ifndef RASTRIGIN_H
#define RASTRIGIN_H

#include "movable.h"
#include <iostream>
using namespace std;

class variable : public abstract_param {
    public:
        static const double VAR_MAX = 5.12;
        static const double VAR_MIN = -5.12;
        double x;
        void init(unsigned int *in_seed);
        void generate_tweak(double theta_bar);
        void restore_tweak();
    private:
        unsigned int *seed;
        double prev_x;
        bool is_restorable;
};

class rastrigin : public movable {
    public:
        rastrigin (int dimension);
        double get_score();
        void print_solution(ostream &o) const;
        ~rastrigin();
    private:
        variable *vars;
        unsigned int seed;

};

#endif
