

#include <iostream>
#include <sstream>
#include <mpi.h>
#include <unistd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "unirandom.h"
#include "plsa.h"
#include "tempCountP.h"
#include "dynDebug.h"
#include "mix/intervalMix.h"
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
    char c;
    bool isprolix = false;
    bool isverbose = false;
    bool issteplog = true;
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
            cerr << "Unrecognized option: " << c << endl;
            return 1;
        }
    }
    char *docname = argv[optind];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;
    unirand48 rnd(mpi.rank);
    fly_params flyParams = readFlyParams(docroot);
    fly theFly(flyParams);
    expHoldP::Param scheParam(docroot);
    tempCountP::Param frozenParam(docroot);
    intervalMix<fly>::Param mixParam(docroot);
    // parallelFBMove<fly, debugSTD, adaptMix> *fly_problem =
    //        new parallelFBMove<fly, debugSTD, adaptMix>(theFly, rnd, docroot, mpi);
    // plsa *pschedule = new plsa(docroot, mpi);
    pannealer<fly, expHoldP, tempCountP, parallelFBMove, intervalMix>
            *fly_sa = new pannealer<fly, expHoldP, tempCountP,
                                    parallelFBMove, intervalMix>
            (theFly, rnd, scheParam, frozenParam, mixParam, docroot, mpi);
    string outprefix = flyParams.infile_name + "_" +
            ((ostringstream*)&(ostringstream()<<mpi.rank))->str();
    if (isprolix) {
        fly_sa->setProlix(file, (outprefix + ".prolix").c_str());
    }

    if (isverbose) {
        fly_sa->setMixLog(file, (outprefix + ".mixlog").c_str());
    }

    if (mpi.rank == 0) {
        fly_sa->setCoolLog(file,(flyParams.infile_name + ".log").c_str());
        // fly_sa->setProlix(file, (flyParams.infile_name + ".prolix").c_str());
    }
    if (issteplog) {
        fly_sa->setStepLog(file, (outprefix + ".steplog").c_str());
    }


    cout << "The initial energy is " << theFly.get_score() << endl;
    fly_sa->loop();
    cout << "The final energy is " << theFly.get_score() << endl;
    if (fly_sa->getWinner() == mpi.rank) {
        theFly.writeAnswer("eqparms");
        fly_sa->writeResult(docroot);
        boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
        write_xml(docname, pt, std::locale(), settings);
    }
    delete fly_sa;
    MPI_Finalize();

    return 0;
}
