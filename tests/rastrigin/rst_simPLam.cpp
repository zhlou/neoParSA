/*
 * rst_simPLam.cpp
 *
 *  Created on: Jul 9, 2014
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
#include "rastrigin.h"

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
    rastrigin rst(xmlroot, rnd);

    simPLam::Param scheParam(xmlroot);
    globalCount::Param frozenParam(xmlroot);
    intervalMix<rastrigin>::Param mixParam(xmlroot);
    pannealer<rastrigin, simPLam, globalCount, parallelFBMove, intervalMix>
            *rst_sa = new pannealer<rastrigin, simPLam, globalCount,
                                    parallelFBMove, intervalMix>
            (rst, rnd, scheParam, frozenParam, mixParam, xmlroot, mpi);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);
    string outprefix = basename+ "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str();
    if (isprolix) {
        rst_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

    if (isverbose) {
        rst_sa->setMixLog(file, (outprefix + ".mixlog").c_str());
    }

    if (mpi.rank == 0) {
        rst_sa->setCoolLog(file,(basename+".log").c_str());
    }

    if (issteplog) {
        rst_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }
    cout << "The initial energy is " << rst.get_score() << endl;
    rst_sa->loop();
    cout << "Final energy is " << rst.get_score() << endl;
    if (rst_sa->getWinner() == mpi.rank) {
        rst.write_section((xmlChar *) "tour");
        rst_sa->writeResult();
        xmlSaveFormatFile(docname, xmldoc,1);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    delete rst_sa;
    MPI_Finalize();

}
