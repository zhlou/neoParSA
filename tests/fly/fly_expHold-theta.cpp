/*
 * fly_expHold-theta.cpp
 *
 *  Created on: Sep 17, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <sstream>
#include <libxml/parser.h>
#include <mpi.h>

#include "pannealer.h"
#include "move/FBNCTheta.h"
#include "unirandom.h"
#include "plsa.h"
#include "tempCountP.h"
#include "dynDebug.h"
#include "adaptMix.h"
#include "intervalMix.h"
#include "pulseNoAdopt.h"
#include "fly.h"
#include "expHoldP.h"

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
    expHoldP::Param scheParam(docroot);
    tempCountP::Param frozenParam(docroot);
    pulseNoAdopt<fly>::Param mixParam(docroot);
    pannealer<fly, expHoldP, tempCountP, FBNCTheta, pulseNoAdopt>
            *fly_sa = new pannealer<fly, expHoldP, tempCountP,
                                    FBNCTheta, pulseNoAdopt>
            (theFly, rnd, scheParam, frozenParam, mixParam, docroot, mpi);
    if (mpi.rank == 0) {
        fly_sa->setCoolLog(file,(flyParams.infile_name + ".log").c_str());
        fly_sa->setProlix(file, (flyParams.infile_name + ".prolix").c_str());
    }
    fly_sa->setStepLog(file, (flyParams.infile_name + "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str() +
            ".steplog").c_str());


    cout << "The initial energy is " << theFly.get_score() << endl;
    fly_sa->loop();
    cout << "The final energy is " << theFly.get_score() << endl;
    if (fly_sa->getWinner() == mpi.rank) {
        theFly.writeAnswer("eqparms");
        fly_sa->writeResult();
        xmlSaveFormatFile(docname, doc, 1);
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();
    delete fly_sa;
    MPI_Finalize();

    return 0;
}







