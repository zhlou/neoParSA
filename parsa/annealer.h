#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
class movable;
/*
 * This should be the abstract annealer class all the actual annealer
 * derived from. It has a constructor with basic xml interpertation function.
 * The simple implemenation with exponential cooling and rejection conting
 * frozen condition is now a separate simpleAnnealer class.
 */
class annealer
{
public:
    annealer(movable *theproblem, xmlNode *root);
    virtual ~annealer();
    double loop();
protected:
    virtual void updateStep(bool accept, double delta) = 0;
    virtual bool frozen() = 0;
    virtual void cool_s() = 0;
    virtual bool inSegment() = 0;
    virtual void updateSegment() = 0;
    xmlNode *xmlroot;
    xmlNode *xmlsection;
    double s;
    double lambda;
    long unsigned step_cnt;
    unsigned rnd_seed;
    movable *problem;
    double energy;
};

#endif
