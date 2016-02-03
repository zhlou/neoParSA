/*
 * rst_pEH-fbMix.cpp
 *
 *  Created on: Aug 27, 2014
 *      Author: zhlou
 */

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <unistd.h>
#include <libgen.h>
#include <getopt.h>
#include <mpi.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#include <boost/optional.hpp>

#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "unirandom.h"
#include "criCountP.h"
#include "dynDebug.h"
#include "mix/feedbackMix.h"
#include "expHoldP.h"

#include "udrst.h"


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
    int iscoollog = 0;
    int readInitStates = 0;
    int optIndex;
    char *stateListFile = NULL;
    struct option long_options[] = {
        {"read-state", 1, &readInitStates, 1},
        {"cool-log", 0, &iscoollog, 1},
        {0, 0, 0, 0}
    };

    std::string section;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ( (c = getopt_long(argc, argv, "pvL", long_options, &optIndex)) != -1) {
            switch(c) {
            case 0:
                switch (optIndex) {
                case 0:
                    stateListFile = optarg;
                    break;
                case 1:
                    break;
                default:
                    throw std::runtime_error("Unrecognized option");
                }
                break;
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
    std::string docname(argv[optind]);
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirandom rnd(mpi.rank);
    boost::optional<unsigned int> seed = docroot.get_optional<unsigned int>("<xmlattr>.seed");
    if (seed) {
        rnd.setSeed(*seed + mpi.rank);
    }
    udrst rst(docroot, rnd);
    expHoldP::Param scheParam(docroot);
    criCountP::Param frozenParam(docroot);
    feedbackMix<udrst>::Param mixParam(docroot);
    pannealer<udrst, expHoldP, criCountP, parallelFBMove, feedbackMix>
            *rst_sa = new pannealer<udrst, expHoldP, criCountP,
                                    parallelFBMove, feedbackMix>
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
        if (iscoollog)
            rst_sa->setCoolLog(file,(basename + ".log").c_str());
        if (isverbose) {
            rst_sa->setMixLog(file, (basename + ".mixlog").c_str());
        }
        // fly_sa->setProlix(file, (flyParams.infile_name + ".prolix").c_str());
    }
    if (issteplog) {
        rst_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }
    if (readInitStates) {
        std::string readStatePrefix;
        std::string line;
        std::ifstream is(stateListFile);
        int i = 0;
        while (!(std::getline(is,line)).eof()) {
            if (mpi.rank == i) {
                readStatePrefix = line;
                break;
            }
            ++i;
        }

        if (!readStatePrefix.empty()) {
            rst_sa->readUnifiedInitState(readStatePrefix);
        } else {
            throw std::runtime_error("unable to find state");
        }
        is.close();
    }

    std::cout << "The initial energy is " << rst.get_score() << std::endl;
    rst_sa->loop();
    std::cout << "The final energy is " << rst.get_score() << std::endl;
    if (rst_sa->getWinner() == mpi.rank) {
        rst.write_section(docroot, "output");
        rst_sa->writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }

    delete rst_sa;
    MPI_Finalize();
    return 0;


}
