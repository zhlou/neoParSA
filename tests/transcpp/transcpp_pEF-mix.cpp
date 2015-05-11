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

#include <unistd.h>
#include <getopt.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <mpi.h>
#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "intervalMix.h"
#include "unirandom.h"
#include "expHoldP.h"
#include "criCountP.h"
#include "dynDebug.h"

using namespace std;

using boost::property_tree::ptree;


int mode_verbose;

int main(int argc, char** argv)
{
    MPIState mpiState;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiState.nnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiState.rank);
    MPI_Comm_group(MPI_COMM_WORLD, &mpiState.group);
    mpiState.comm = MPI_COMM_WORLD;

    bool isprolix = false;
    bool issteplog = false;
    bool isverbose = false;
    int iscoollog = 0;
    int readInitStates = 0;
    int optIndex;
    char *stateListFile = NULL;
    const char *readStatePrefix = NULL;
    struct option long_options[] = {
        {"read-state", 1, &readInitStates, 1},
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
    rnd.setSeed(seed+mpiState.rank);
    //rnd.setSeed(getpid());

    xmlDoc *doc = xmlParseFile(xmlname.c_str());
    xmlNode *docroot = xmlDocGetRootElement(doc);
    expHoldP::Param scheParam(docroot);
    criCountP::Param frozenParam(docroot);
    intervalMix<Organism>::Param mixParam(docroot);
    pannealer<Organism, expHoldP, criCountP, parallelFBMove, intervalMix>
            *annealer = new pannealer<Organism, expHoldP, criCountP,
            parallelFBMove, intervalMix>(embryo, rnd, scheParam, frozenParam,
            mixParam, docroot, mpiState);
    string bname(xmlname);
    size_t sz=bname.size();
    bname.resize(sz-4);
    string outprefix = bname + "_" +
            ((ostringstream*)&(ostringstream() << mpiState.rank))->str();

    if (0 == mpiState.rank) {
        if (iscoollog)
            annealer->setCoolLog(file, (bname + ".log").c_str());
        if (isprolix)
            annealer->setProlix(file, (bname + ".prolix").c_str());
        if (isverbose) {
            annealer->setMixLog(file, (bname + ".mixlog").c_str());
        }
    }

    if (issteplog)
        annealer->setStepLog(file, (outprefix + ".steplog").c_str());


    if (readInitStates) {
        std::string line;
        std::ifstream is(stateListFile);
        int i = 0;
        while (!(std::getline(is,line)).eof()) {
            if (mpiState.rank == i) {
                readStatePrefix = line.c_str();
                break;
            }
            ++i;
        }
        if (readStatePrefix) {
            annealer->readUnifiedInitState(readStatePrefix);
        } else {
            throw std::runtime_error("unable to find state");
        }
        is.close();
    }

    cerr << "The energy is " << embryo.get_score() << endl;
    annealer->loop();
    cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
    xmlFreeDoc(doc);

    //embryo.printParameters(cerr);

    if (annealer->getWinner() == mpiState.rank) {
        embryo.write("Output", root_node);
        annealer->ptreeGetResult(root_node);
        boost::property_tree::xml_writer_settings<char> settings(' ', 2);
        write_xml(xmlname, pt, std::locale(), settings);
    }
    xmlCleanupParser();
    delete annealer;
    MPI_Finalize();


    return 0;
}

