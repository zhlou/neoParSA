#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
#include "unirandom.h"
#include "aState.h"
#include "dynDebug.h"
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
template <class Problem, class Schedule, class FrozenCnd, template<class> class Move>
class annealer
{
public:
    annealer(Problem &problem, unirandom * const in_rand,
             typename Schedule::Param scheParam,
             typename FrozenCnd::Param frozenParam, xmlNode *root);
    virtual ~annealer();
    double loop();
    double initMoves();
    // fixed T moves with no stats updated for schedule. move part is involved
    // hence likely updated.
    double fixedTMoves(double S, long steps);
    void writeResult();
    void setStepLog(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st, outname); debugOut.precision(16);}
    void setCoolLog(debugStatus st, const char*outname=NULL)
    {cooling->setDebug(st, outname);}
    void setProlix(debugStatus st, const char*outname=NULL)
    {move->setDebug(st, outname);}
protected:
    dynDebug debugOut;
    Schedule *cooling;
    FrozenCnd *frozen;
    Move<Problem> *move;
    unirandom *const rand;
    xmlNode *xmlroot;
    aState state;
    int initLoop;
    double initS;
    bool is_init;
    double tlaps;
    annealer(unirandom * const in_rand, xmlNode *root);
    virtual void updateSegment(aState &state) {cooling->updateSegment(state);}

    bool step();
    void initState(xmlNode* root);
    virtual void writeResultData(xmlNode* result);
    virtual void writeMethodText(xmlNode* method);
};


#include "annealer.hpp"

#endif
