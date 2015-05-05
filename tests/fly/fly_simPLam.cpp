/*
 * fly_simPLam.cpp
 *
 *  Created on: Dec 18, 2013
 *      Author: zhlou
 */




#include <iostream>
#include <sstream>
#include <libxml/parser.h>
#include <unistd.h>
#include <mpi.h>

#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "unirandom.h"
#include "simPLam.h"
#include "globalCount.h"
#include "dynDebug.h"
#include "intervalMix.h"
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
    char c;
    bool isprolix = false;
    bool isverbose = false;
    bool issteplog = false;

    while ( (c = getopt(argc, argv, "lpv")) != -1) {
        switch(c) {
        case 'l':
            issteplog = true;
            break;
        case 'p':
            isprolix = true;
            break;
        case 'v':
            isverbose = true;
            break;
        default:
            cerr << "Unrecognized option: " << c << endl;
            return 1;
        }
    }
    char *docname = argv[optind];
    xmlDoc *doc;
    if ( (doc = xmlParseFile(docname) ) == NULL) {
        cerr << "XML file " << docname << " ill formed or not found!" << endl;
        return 2;
    }
    xmlNode *docroot;
    if ( (docroot = xmlDocGetRootElement(doc) ) == NULL) {
        cerr << "Input incorrect" << endl;
        return 2;
    }
    unirand48 rnd(mpi.rank);
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    simPLam::Param schdParam(docroot);
    globalCount::Param frozenParam(docroot);
    intervalMix<fly>::Param mixParam(docroot);

    pannealer<fly, simPLam, globalCount, parallelFBMove, intervalMix>
        *fly_sa = new pannealer<fly, simPLam, globalCount, parallelFBMove,
                                intervalMix>
                    (theFly, rnd, schdParam, frozenParam, mixParam, docroot, mpi);
    string outprefix = flyParams.infile_name + "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str();
    if (isprolix) {
        fly_sa->setProlix(file, (outprefix+".prolix").c_str());
    }
    if (isverbose) {
        fly_sa->setMixLog(file, (outprefix + ".mixlog").c_str());
    }
    if (mpi.rank == 0) {
        fly_sa->setCoolLog(file, (flyParams.infile_name + ".log").c_str());
    }
    if (issteplog) {
        fly_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }

    fly_sa->loop();
    if (mpi.rank == (fly_sa->getWinner())) {
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
