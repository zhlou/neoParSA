#ifndef ANNEALER_H
#define ANNEALER_H

#include "movable.h"
#include <libxml/parser.h>

class annealer {
    public:
        annealer(movable *theproblem, xmlNode *root);
        ~annealer();
        double loop();
    private:
        xmlNode *xmlroot;
        double s;
        bool frozen();
        unsigned reject_cnt;
        unsigned rnd_seed;
        movable *problem;
};


#endif
