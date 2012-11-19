#ifndef ANNEALER_H
#define ANNEALER_H

#include "movable.h"

class annealer {
    public:
        annealer(movable *theproblem);
        ~annealer();
        double loop();
    private:
        double s;
        bool frozen();
        unsigned reject_cnt;
        unsigned rnd_seed;
        movable *problem;
};


#endif
