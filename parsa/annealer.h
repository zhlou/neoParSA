#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
class movable;
class annealer {
    public:
        annealer(movable *theproblem, xmlNode *root);
        virtual ~annealer();
        double loop();
    protected:
        xmlNode *xmlroot;
        double s;
        double lambda;
        virtual bool frozen();
        virtual void cool_s();
        unsigned reject_cnt;
        unsigned rnd_seed;
        movable *problem;
};


#endif
