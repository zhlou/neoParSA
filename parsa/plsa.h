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
#include "annealer.h"

class plsa: public lam
{
public:
    plsa(movable *theproblem, xmlNode *root, MPI_Comm thecomm, int in_nnodes,
            int in_rank);
    ~plsa();
protected:
    struct StatData {
        double s;
        double mean;
        double sd;
        double energy;
        long success;
        //long moves;
    };
    //void updateS();
    void initStats();
    void updateSegment();
    bool frozen();
    void CommSegment();
    bool CommCheckFrozen();
    StatData l_stat, *local_stat_buf;
    MPI_Win stat_win;
    MPI_Group group;
    MPI_Comm comm;

    int nnodes;
    int rank;
};

#endif /* PLSA_H_ */
