#include <iostream>
#include <unistd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#include <boost/optional.hpp>

#include "DoS/DoS.h"
#include "DoS/PWLE.h"
#include "move/feedbackMove.h"
#include "unirandom.h"

#include "fly.h"


int main(int argc, char ** argv)
{
    MPIState mpi;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpi.group);
    mpi.comm = MPI_COMM_WORLD;

    char *outname = NULL;
    char *readfile = NULL;
    char *savefile = NULL;
    char c;
    while ( (c = getopt(argc, argv, "o:s:r:")) != -1) {
        switch(c) {
        case 'o':
            outname = optarg;
            break;
        case 'r':
            readfile = optarg;
            break;
        case 's':
            savefile = optarg;
            break;
        default:
            std::cerr << "Unrecognized option " << c << std::endl;
            break;
        }
    }
    if (argc <= optind) {
        std::cerr << "Missing input files" << std::endl;
        std::cerr << "Usage: fly_wls [-o output] [-s savefile] [-r readfile] input " << std::endl;
        return 1;
    }
    char *docname = argv[optind];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &xmlroot = pt.begin()->second;

    unirand48 rnd(mpi.rank);
    fly_params flyParams = readFlyParams(xmlroot);
    fly theFly(flyParams);
    feedbackMove<fly> flyMove(theFly, rnd, xmlroot);
    DoS<fly, feedbackMove, PWLE>::Param param;
    param.estParam.mpi = &mpi;

    ptree &dos_attr = xmlroot.get_child("DoS.<xmlattr>");
    param.initWeight = dos_attr.get<double>("weight", 1e-2);
    param.nSteps = dos_attr.get<long>("nsteps", 100000);

    ptree &wle_attr = xmlroot.get_child("WLEstimator.<xmlattr>");
    param.estParam.eMin = wle_attr.get<double>("eMin", 0);
    param.estParam.binWidth = wle_attr.get<double>("binWidth", 100);
    param.estParam.nBins = wle_attr.get<unsigned int>("nBins", 100000);
    param.estParam.syncFreq = wle_attr.get<unsigned int>("syncFreq", 1000);
    param.estParam.saveFreq = wle_attr.get<unsigned int>("saveFreq", 0);
    param.estParam.saveName = wle_attr.get<std::string>("saveName").c_str();
    
    DoS<fly, feedbackMove, PWLE> simulate(theFly, flyMove, rnd, param);
    PWLE &estm = simulate.getEstimator();
    if (readfile) {
        estm.readHist(readfile);
    }
    simulate.estimate();
    estm.syncHist();
    if (0 == mpi.rank) {
        if (outname == NULL)
            estm.printHist(std::cout);
        else {
            ofstream myout(outname);
            estm.printHist(myout);
            myout.close();
        }
        if (savefile) {
            estm.saveHist(savefile);
        }
    }
    MPI_Finalize();
    return 0;


}
