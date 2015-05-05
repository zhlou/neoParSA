#include "flags.h"
#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "conc.h"
#include "twobit.h"
#include "organism.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <unistd.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <mpi.h>
#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "intervalMix.h"
#include "unirandom.h"
#include "expHoldP.h"
#include "tempCountP.h"
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
    bool iscoollog = true;
    std::string binname(basename(argv[0]));
    try {
        char c;
        while ((c = getopt(argc, argv, "plvN")) != -1) {
            switch (c) {
            case 'l':
                issteplog = true;
                break;
            case 'N':
                iscoollog = false;
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

    ptree & input = pt.get_child("Input");

    Organism embryo(input);
    //embryo.printParameters(cerr);

    unirand48 rnd(mpiState.rank);
    //rnd.setSeed(getpid());

    xmlDoc *doc = xmlParseFile(xmlname.c_str());
    xmlNode *docroot = xmlDocGetRootElement(doc);
    expHoldP::Param scheParam(docroot);
    tempCountP::Param frozenParam(docroot);
    intervalMix<Organism>::Param mixParam(docroot);
    pannealer<Organism, expHoldP, tempCountP, parallelFBMove, intervalMix>
            *annealer = new pannealer<Organism, expHoldP, tempCountP,
            parallelFBMove, intervalMix>(embryo, rnd, scheParam, frozenParam,
            mixParam, docroot, mpiState);
    string outprefix = xmlname + "_" +
            ((ostringstream*)&(ostringstream() << mpiState.rank))->str();

    if (0 == mpiState.rank) {
        if (iscoollog)
            annealer->setCoolLog(file, (xmlname + ".log").c_str());
        if (isprolix)
            annealer->setProlix(file, (xmlname + ".prolix").c_str());
        if (isverbose) {
            annealer->setMixLog(file, (xmlname + ".mixlog").c_str());
        }
    }

    if (issteplog)
        annealer->setStepLog(file, (outprefix + ".steplog").c_str());



    cerr << "The energy is " << embryo.get_score() << endl;
    annealer->loop();
    cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
    xmlFreeDoc(doc);

    //embryo.printParameters(cerr);

    if (annealer->getWinner() == mpiState.rank) {
        ptree output, anneal_output;
        embryo.write("Output", output);
        annealer->ptreeGetResult(output.get_child("Output"));
        //ptree& opt = output.get_child("Output");
        output.put_child("Output.anneal_output", anneal_output);
        boost::property_tree::xml_writer_settings<char> settings(' ', 2);
        write_xml_element(infile, basic_string<ptree::key_type::value_type>(), output, -1, settings);
    }
    xmlCleanupParser();
    delete annealer;
    MPI_Finalize();


    return 0;
}

