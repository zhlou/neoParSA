#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "datatable.h"
#include "twobit.h"
#include "organism.h"
#include "mode.h"
#include "utils.h"
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

#include <iostream>
#include <unistd.h>
#include "DoS/DoS.h"
#include "DoS/PWLE.h"
#include "move/feedbackMove.h"
#include "unirandom.h"


using namespace std;

using boost::property_tree::ptree;


int mode_verbose;



int main(int argc, char **argv)
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
        std::cerr << "Usage: udrstWLS [-o output] [-s savefile] [-r readfile] input " << std::endl;
        return 1;
    }
    mode_verbose = 0;
    string xmlname(argv[optind]);
    fstream infile(xmlname.c_str());

    ptree pt;
    read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);

    ptree& root_node = pt.get_child("Root");
    ptree& mode_node = root_node.get_child("Mode");
    ptree& input_node = root_node.get_child("Input");
    mode_ptr mode(new Mode(xmlname, mode_node));

    Organism embryo(input_node, mode);
    //embryo.printParameters(cerr);
    unsigned int seed = mode->getSeed();
    if (mode->getVerbose() >= 1)
        cerr << "Beginning annealing with seed " << seed << endl;
    unirand48 rnd;
    rnd.setSeed(seed+mpi.rank);

    feedbackMove<Organism> transcMove(embryo, rnd, root_node);
    DoS<Organism, feedbackMove, PWLE>::Param param;

    param.initWeight = 1e-2;
    param.nSteps = 100000;
    param.estParam.eMin=0;
    param.estParam.binWidth=0.01;
    param.estParam.nBins=4040;
    param.estParam.syncFreq = 1000;
    param.estParam.saveName = NULL;
    param.estParam.saveFreq = 0;
    param.estParam.mpi = &mpi;
    boost::optional<ptree &> dos_attr = root_node.get_child_optional("DoS.<xmlattr>");
    if (dos_attr) {
        param.initWeight = dos_attr->get<double>("weight", 1e-2);
        param.nSteps = dos_attr->get<long>("nsteps", 100000);
    }

    boost::optional<ptree &> wle_attr = root_node.get_child("WLEstimator.<xmlattr>");
    {
        param.estParam.eMin = wle_attr->get<double>("eMin", 0);
        param.estParam.binWidth = wle_attr->get<double>("binWidth", 0.01);
        param.estParam.nBins = wle_attr->get<unsigned int>("nBins", 4040);
        param.estParam.syncFreq = wle_attr->get<unsigned int>("syncFreq", 1000);
        param.estParam.saveFreq = wle_attr->get<unsigned int>("saveFreq", 0);
        param.estParam.saveName = wle_attr->get<std::string>("saveName").c_str();
    }

    DoS<Organism, feedbackMove, PWLE> simulate(embryo, transcMove, rnd, param);
    PWLE &estm=simulate.getEstimator();
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
