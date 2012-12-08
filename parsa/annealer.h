#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
class movable;
class annealer {
    public:
        annealer(movable *theproblem, xmlNode *root);
        virtual ~annealer();
        virtual double loop();
    protected:
        xmlNode *xmlroot;
        xmlNode *xmlsection;
        double s;
        double lambda;
        double move();
        unsigned reject_cnt;
        unsigned rnd_seed;
        movable *problem;
};


#endif
