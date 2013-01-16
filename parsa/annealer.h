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
template <class Schedule, class Move>
class annealer
{
public:
    annealer(movable *theproblem, xmlNode *root);
    virtual ~annealer();
    double loop();
    double initMoves();
protected:
    Schedule cooling;
    Move control;
    virtual void updateInitStep(bool accept);
    virtual void updateStep(bool accept) = 0;
    virtual void initStats() = 0;
    virtual bool frozen() = 0;
    virtual void updateS() = 0;
    virtual bool inSegment() = 0;
    virtual void updateSegment() = 0;
    virtual void resetSegmentStats();
    bool move();
    xmlNode *xmlroot;
    xmlNode *xmlsection;
    double s;
    double lambda;
    long unsigned step_cnt;
    unsigned rnd_seed;
    movable *problem;
    double energy;
    int init_loop;
    bool is_init;
};

#include "annealer.hpp"

#endif
