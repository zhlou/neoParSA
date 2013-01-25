#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
//class movable;
/*
 * This is the abstract annealer template. All the actual functionalities
 * come from Schedule and Move class.
 *
 * Schedule class is the (stateful) class that controls the cooling schedule
 * and frozen condition. It needs to have the following interfaces
 *
 * class Schedule
 * {
 * public:
 *      double getInitS();
 *      int getInitLoop();
 *      void updateInitStep(bool accept, double energy);
 *      void initStats();
 *      bool frozen();
 *      void resetSegmentStats();
 *      void updateStep(bool accept, double energy);
 *      double updateS(double s); //returns the updated s
 *      bool inSegment();
 *      void updateSegment();
 * };
 *
 * Move should be a template with move generation and move control that has
 * the actual problem as the parameter. It needs to have the following
 * interfaces
 *
 * template <class Problem>
 * class Move
 * {
 * public:
 *      double get_score();
 *      double propose(); // returns the difference of energies
 *      void accept();
 *      void reject();
 * };
 */
template <class Schedule, class Move, class Random>
class annealer
{
public:
    annealer(Schedule &in_cool, Move &in_move, Random &in_rand, xmlNode *root);
    virtual ~annealer();
    double loop();
    double initMoves();
protected:
    Schedule &cooling;
    Move &move;
    Random &rand;
    //virtual void updateInitStep(bool accept);
    //virtual void updateStep(bool accept) = 0;
    //virtual void initStats() = 0;
    //virtual bool frozen() = 0;
    //virtual void updateS() = 0;
    //virtual bool inSegment() = 0;
    //virtual void updateSegment() = 0;
    //virtual void resetSegmentStats();
    bool step();
    xmlNode *xmlroot;
    //xmlNode *xmlsection;
    double s;
    //double lambda;
    long unsigned step_cnt;
    int initLoop;
    double initS;

    //unsigned rnd_seed;
    //movable *problem;
    double energy;
    //int init_loop;
    bool is_init;
};



#include "annealer.hpp"

#endif
