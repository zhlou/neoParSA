/*
 * rst_pEH-best.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <unistd.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <mpi.h>

#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "unirandom.h"
#include "tempCountP.h"
#include "dynDebug.h"
#include "periodBest.h"
#include "expHoldP.h"

#include "rastrigin.h"


int main(int argc, char **argv)
{
    MPIState mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpi.group);
    mpi.comm = MPI_COMM_WORLD;
    bool isprolix = false;
    bool isverbose = false;
    bool issteplog = true;

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt(argc, argv, "pvL")) != -1) {
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
        std::cerr << "Usage: " << binname << " [ -x section_name ] input_file"
                  << std::endl;
        return -1;
    }
    char *docname = argv[optind];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (docroot == NULL) {
        std::cerr << "Input incorrect" << std::endl;
        return -1;
    }
    unirandom rnd(mpi.rank);
    rastrigin rst(docroot, rnd);
    expHoldP::Param scheParam(docroot);
    tempCountP::Param frozenParam(docroot);
    periodBest<rastrigin>::Param mixParam(docroot);
    pannealer<rastrigin, expHoldP, tempCountP, parallelFBMove, periodBest>
            *rst_sa = new pannealer<rastrigin, expHoldP, tempCountP,
                                    parallelFBMove, periodBest>
            (rst, rnd, scheParam, frozenParam, mixParam, docroot, mpi);
    string basename(docname);
    size_t sz = basename.size();
    basename.resize(sz-4);

    string outprefix = basename + "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str();
    if (isprolix) {
        rst_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }



    if (mpi.rank == 0) {
        rst_sa->setCoolLog(file,(basename + ".log").c_str());
        if (isverbose) {
            rst_sa->setMixLog(file, (basename + ".mixlog").c_str());
        }
        // fly_sa->setProlix(file, (flyParams.infile_name + ".prolix").c_str());
    }
    if (issteplog) {
        rst_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }

    std::cout << "The initial energy is " << rst.get_score() << std::endl;
    rst_sa->loop();
    std::cout << "The final energy is " << rst.get_score() << std::endl;
    if (rst_sa->getWinner() == mpi.rank) {
        rst.write_section((xmlChar *)"output");
        rst_sa->writeResult();
        xmlSaveFormatFile(docname, doc, 1);
    }

    xmlFreeDoc(doc);
    xmlCleanupParser();
    delete rst_sa;
    MPI_Finalize();
    return 0;


}


