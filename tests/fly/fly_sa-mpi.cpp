/*
 * fly_sa-mpi.cpp
 *
 *  Created on: Mar 3, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <libxml/parser.h>
#include <mpi.h>

#include "pannealer.h"
#include "parallelFBMove.h"
#include "unirandom.h"
#include "plsa.h"
#include "dynDebug.h"
#include "adaptMix.h"
#include "fly.h"

using namespace std;

int main(int argc, char **argv)
{
    MPIState mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpi.group);
    mpi.comm = MPI_COMM_WORLD;

    if (argc <= 1) {
        cerr << "Missing input files" << endl;
        return 1;
    }
    char *docname = argv[1];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (docroot == NULL) {
        cerr << "Input incorrect" << endl;
        return 2;
    }
    unirand48 rnd(mpi.rank);
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    // parallelFBMove<fly, debugSTD, adaptMix> *fly_problem =
    //        new parallelFBMove<fly, debugSTD, adaptMix>(theFly, rnd, docroot, mpi);
    // plsa *pschedule = new plsa(docroot, mpi);
    pannealer<fly, plsa, parallelFBMove, adaptMix> *fly_sa =
            new pannealer<fly, plsa, parallelFBMove, adaptMix>(theFly, &rnd, docroot, mpi);
    cout << "The initial energy is " << theFly.get_score() << endl;
    fly_sa->loop();
    cout << "The final energy is " << theFly.get_score() << endl;
    if (fly_sa->getWinner() == mpi.rank)
        theFly.writeAnswer();
    delete fly_sa;
    MPI_Finalize();

    return 0;
}




