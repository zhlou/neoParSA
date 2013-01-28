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

class plsa: public lam
{
public:
    plsa(xmlNode *root, MPI_Comm thecomm, int in_nnodes,
            int in_rank);
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
    MPI_Group group;
    MPI_Comm comm;

    int nnodes;
    int rank;
};

#endif /* PLSA_H_ */
