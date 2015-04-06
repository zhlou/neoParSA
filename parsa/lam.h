/*
 * lam.h
 *
 *  Created on: Dec 8, 2012
 *      Author: zhlou
 */

#ifndef LAM_H_
#define LAM_H_

//#include "annealer.h"

#include <libxml/tree.h>
#include "aState.h"
#include "invLinearFit.h"
#include "dynDebug.h"

class lam
{
public:
    class Param {
    public:
        int proc_tau;
        //double freeze_crit;
        //int cnt_crit;
        double lambda;
        double w_mean;
        double w_sd;
        debugStatus st;
        const char * outname;
        Param(xmlNode *root, debugStatus st=ignore, const char *outname=NULL);
    };
    lam(Param param);
    // lam(xmlNode *root);
    virtual ~lam();
    //double getInitS();
    //int getInitLoop();
    //bool frozen(const aState state);
    void updateInitStep(bool accept, aState state);
    void initStats(aState state);
    void initStats(double initMean, double initVar, double initAccRatio, aState state);

    void updateStep(bool accept, aState state);
    double updateS(aState state);
    bool inSegment(aState state){return !((state.step_cnt % proc_tau) == 0);}
    void updateStats(aState state);

//    virtual bool needMix(){return false;}
    void setDebug(debugStatus st, const char*outname=NULL)
        {debugOut.setDebug(st, outname);}
    static const char* name;

protected:
    // parameters
    int proc_tau;
    //double freeze_crit;
    //int cnt_crit;
    dynDebug debugOut;
    double lambda;
    double w_mean;
    double w_sd;


    // states
    invLinearFit *fit_mean, *fit_sd;
    double alpha;
    double acc_ratio;
    double vari;
    double mean;
    int success;
    double old_energy;
    int freeze_cnt;
    int tau_count;

    //double energy;
    //double s;



    static const double UNINITIALIZED;

    void updateSegment(aState state);
    void resetSegmentStats();
    
    void collectStats();
    virtual void collectInitStats(unsigned long init_loop);
    virtual void collectInitStats(double initMean, double initVar, double initAccRatio);
    void resetLam();
    void updateLam();
    void initStatsCore(const aState &state);
    virtual void updateEstimators(double s);
    //void local_frozen(const aState& state);
    //virtual bool global_frozen();



private:

    //int init_loop;
    //double init_S;
    //long step_cnt;

};

#endif /* LAM_H_ */
