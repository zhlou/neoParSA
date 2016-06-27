/*
 * fly_expHold.cpp
 *
 *  Created on: Jun 13, 2013
 *      Author: zhlou
 */

#include <iostream>
#include <unistd.h>
#include <getopt.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "pannealer.h"
#include "move/parallelCauchyMove.h"
#include "mix/feedbackMix.h"
#include "unirandom.h"
#include "expHoldP.h"
#include "criCountP.h"
#include "fly.h"

int main(int argc, char **argv)
{
    MPIState mpiState;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiState.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiState.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpiState.group);
    mpiState.comm = MPI_COMM_WORLD;

    bool isprolix = false;
    bool issteplog = false;
    bool ismixlog = false;
    int iscoollog = 0;
    int readInitState = 0;
    char *stateListFile = NULL;
    int optIndex;

    struct option long_options[] = {
        {"read-state", 1, &readInitState, 1},
        {"cool-log", 0, &iscoollog, 1},
        {0, 0, 0, 0}
    };

    std::string binname(basename(argv[0]));
    try {
        char c;
        while ((c = getopt_long(argc, argv, "plv", long_options, &optIndex)) != -1) {
            switch (c) {
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
            case 'l':
                issteplog = true;
                break;
            case 'p':
                isprolix = true;
                break;
            case 'v':
                ismixlog = true;
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
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirand48 rnd(mpiState.rank);
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);

    expHoldP::Param scheduleParam(docroot);
    criCountP::Param frozenParam(docroot);
    feedbackMix<fly>::Param mixParam(docroot);
    pannealer<fly, expHoldP, criCountP, parallelCauchyMove, feedbackMix>
        *sa = new pannealer<fly, expHoldP, criCountP, parallelCauchyMove, feedbackMix>(theFly, rnd, scheduleParam, frozenParam, mixParam, docroot, mpiState);

    if (0 == mpiState.rank) {
        if (isprolix) {
            sa->setProlix(file,(flyParams.infile_name+".prolix").c_str());
        }
        if (iscoollog)
            sa->setCoolLog(file,(flyParams.infile_name+".log").c_str());
        if (ismixlog) {
            sa->setMixLog(file, (flyParams.infile_name+".mixlog").c_str());
        }

    }
    if (issteplog)
        sa->setStepLog(file, (flyParams.infile_name+".steplog").c_str());


    if (readInitState) {
        std::string readStatePrefix;
        std::string line;
        std::ifstream is(stateListFile);
        int i = 0;
        while (!(std::getline(is,line)).eof()) {
            if (mpiState.rank == i) {
                readStatePrefix = line;
                break;
            }
            ++i;
        }
        if (! readStatePrefix.empty()) {
            sa->readUnifiedInitState(readStatePrefix);
        } else {
            throw std::runtime_error("unable to find state");
        }
        is.close();
    }
    sa->loop();
    if (sa->getWinner() == mpiState.rank) {
        theFly.writeAnswer("eqparms");
        sa->writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }
    delete sa;
    MPI_Finalize();

    return 0;
}
