/*
 * tsp_simPLam.cpp
 *
 *  Created on: Jul 3, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <cstring>
#include <libxml/parser.h>
#include <libgen.h>

#include "unirandom.h"
#include "simPLam.h"
#include "globalCount.h"
#include "move/parallelFBMove.h"
#include "pannealer.h"
#include "intervalMix.h"
#include "tsp.h"

using namespace std;

int main(int argc, char **argv)
{
    MPIState mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpi.group);
    mpi.comm = MPI_COMM_WORLD;

    bool isprolix = false;
    bool issteplog = true;
    bool isverbose = false;

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt(argc, argv, "vpL")) != -1) {
            switch(c) {
            case 'L':
                issteplog = false;
                break;
            case 'p':
                isprolix = true;
                break;
            case 'v':
                isverbose = true;
                break;
            default:
                throw std::runtime_error("Unrecognized option");
            }
        }
        if (argc <= optind) {
            throw std::runtime_error("Missing input file");
        }
    } catch (std::exception &ex) {
        std::cerr << ex.what() << std::endl;
        std::cerr << "Usage: " << binname << " [-v] [-L] [-p] input_file"
                  << std::endl;
        return -1;
    }
    char *docname = argv[optind];
    xmlDocPtr xmldoc = xmlParseFile(docname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);

    unirand48 rnd;
    tsp theTSP(xmlroot);

    simPLam::Param scheParam(xmlroot);
    globalCount::Param frozenParam(xmlroot);
    intervalMix<tsp>::Param mixParam(xmlroot);
    pannealer<tsp, simPLam, globalCount, parallelFBMove, intervalMix>
            *tsp_sa = new pannealer<tsp, simPLam, globalCount,
                                    parallelFBMove, intervalMix>
            (theTSP, rnd, scheParam, frozenParam, mixParam, xmlroot, mpi);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);
    string outprefix = basename+ "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str();
    if (isprolix) {
        tsp_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

    if (isverbose) {
        tsp_sa->setMixLog(file, (outprefix + ".mixlog").c_str());
    }

    if (mpi.rank == 0) {
        tsp_sa->setCoolLog(file,(basename+".log").c_str());
    }

    if (issteplog) {
        tsp_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }
    cout << "The initial energy is " << theTSP.get_score() << endl;
    tsp_sa->loop();
    cout << "Final energy is " << theTSP.get_score() << endl;
    cout << "Actual energy is " << theTSP.calc_tour() << endl;
    if (tsp_sa->getWinner() == mpi.rank) {
        theTSP.write_tour(xmlroot, "tour");
        tsp_sa->writeResult();
        xmlSaveFormatFile(docname, xmldoc,1);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    delete tsp_sa;
    MPI_Finalize();

}



