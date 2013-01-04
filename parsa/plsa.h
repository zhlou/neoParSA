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
#include "PStats.h"

class plsa: public lam
{
public:
    plsa(movable *theproblem, xmlNode *root, MPI_Comm thecomm, int in_nnodes,
            int in_rank);
    ~plsa();
protected:
    //void updateS();
    void initStats();
    void updateSegment();
    bool frozen();
    PStats StatsComm;
    MPI_Comm comm;
    int nnodes;
    int rank;
};

#endif /* PLSA_H_ */
