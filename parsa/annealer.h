#ifndef ANNEALER_H
#define ANNEALER_H

#include <libxml/parser.h>
#ifdef USE_BOOST
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#endif
#include "unirandom.h"
#include "aState.h"
#include "dynDebug.h"
#include "onePassMeanVar.h"
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
 *      //double getInitS(); // no longer used
 *      //int getInitLoop(); // no longer used
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
    annealer(Problem &problem, unirandom& in_rand,
             typename Schedule::Param scheParam,
             typename FrozenCnd::Param frozenParam, xmlNode *root);
    virtual ~annealer();
    double loop();
    double initMoves();
    double initMovesOnly();
    // fixed T moves with no stats updated for schedule. move part is involved
    // hence likely updated.
    double fixedTMoves(double S, long steps);
    void saveUnifiedInitState(const char * filename);
    double readUnifiedInitState(const char * filename);
#ifdef USE_BOOST
    virtual void ptreeGetResult(ptree &pt);
#endif
    void writeResult();
    void setStepLog(debugStatus st, const char* outname=NULL)
    {debugOut.setDebug(st, outname); debugOut.precision(16);}
    void setCoolLog(debugStatus st, const char*outname=NULL)
    {cooling->setDebug(st, outname);}
    void setProlix(debugStatus st, const char*outname=NULL)
    {move->setDebug(st, outname);}
protected:
    Problem &problem;
    dynDebug debugOut;
    Schedule *cooling;
    FrozenCnd *frozen;
    Move<Problem> *move;
    unirandom& rand;
    xmlNode *xmlroot;
    aState state;
    int initLoop;
    double initS;
    bool is_init;
    double tlaps;
    
    double initMean;
    double initVar;
    double initAccRatio;
    
    annealer(Problem &problem, unirandom& in_rand, xmlNode *root);
    virtual void updateStats(aState &state) {cooling->updateStats(state);}

    bool step();
    void initState(xmlNode* root);
    virtual void writeResultData(xmlNode* result);
    virtual void writeMethodText(xmlNode* method);
};


#include "annealer.hpp"

#endif
