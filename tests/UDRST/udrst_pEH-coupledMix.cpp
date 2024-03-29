/*
 * File:   udrst_pEH-coupledMix.cpp
 * Author: zhlou
 *
 * Created on February 6, 2015, 5:30 PM
 */

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <string>

#include <unistd.h>
#include <libgen.h>
#include <mpi.h>

#include "coupledAnnealer.h"
#include "FBMoveIntervalMix.h"
#include "unirandom.h"
#include "tempCountP.h"
#include "dynDebug.h"
#include "expHoldP.h"

#include "udrst.h"

using namespace std;

/*
 *
 */
int main(int argc, char** argv)
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
    std::string docname(argv[optind]);
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;

    unirandom rnd(mpi.rank);
    udrst rst(docroot, rnd);
    expHoldP::Param scheParam(docroot);
    tempCountP::Param frozenParam(docroot);
    FBMoveIntervalMix<udrst>::Param mixParam(docroot);
    coupledAnnealer<udrst, expHoldP, tempCountP, FBMoveIntervalMix>
            *rst_sa = new coupledAnnealer<udrst, expHoldP, tempCountP,
                                    FBMoveIntervalMix>
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
        rst.write_section(docroot, "output");
        rst_sa->writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }

    delete rst_sa;
    MPI_Finalize();
    return 0;

}
