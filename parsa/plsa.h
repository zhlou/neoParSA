/*
 * plsa.h
 *
 *  Created on: Dec 6, 2012
 *      Author: zhlou
 */

#ifndef PLSA_H_
#define PLSA_H_
#include <mpi.h>
#include <libxml/tree.h>
#include "lam.h"
#include "MPIState.h"

class plsa: public lam
{
public:
    plsa(xmlNode *root, const MPIState &mpiState);
    ~plsa();
    //bool frozen(aState state);
protected:
    struct StatData {
        //double s;
        double mean;
        double var;
        //double energy;
        long success;
        //long moves;
    };
    //void updateS();
    //void initStats();
    void collectInitStats(unsigned long init_loop);
    void updateEstimators(double s);
    bool global_frozen();

    void PackNCommStats(bool UseSD = true);
//    void initEstimators();

    StatData l_stat, *local_stat_buf;
    MPI_Win stat_win;

    const MPIState &mpi;

};

#endif /* PLSA_H_ */
