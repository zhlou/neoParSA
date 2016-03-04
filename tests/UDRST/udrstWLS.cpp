/*
 * rstWLS.cpp
 */
#include <iostream>
#include <unistd.h>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#include <boost/optional.hpp>

#include "DoS/DoS.h"
#include "DoS/WLEstimator.h"
#include "move/feedbackMove.h"
#include "udrst.h"
#include "unirandom.h"




int main(int argc, char **argv)
{

    unirandom rnd;
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

    fstream infile(argv[optind]);
    ptree pt;
    read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
    infile.close();
    ptree &xmlroot = pt.begin()->second;

    udrst rst(xmlroot, rnd);
    feedbackMove<udrst> rstMove(rst, rnd, xmlroot);
    DoS<udrst, feedbackMove, WLEstimator>::Param param;
    param.initWeight = 1e-2;
    param.nSteps = 100000;
    param.estParam.eMin=0;
    param.estParam.binWidth=0.01;
    param.estParam.nBins=4040;
    boost::optional<ptree &> dos_attr = xmlroot.get_child_optional("DoS.<xmlattr>");
    if (dos_attr) {
        param.initWeight = dos_attr->get<double>("weight", 1e-2);
        param.nSteps = dos_attr->get<long>("nsteps", 100000);
    }

    boost::optional<ptree &> wle_attr = xmlroot.get_child("WLEstimator.<xmlattr>");
    {
        param.estParam.eMin = wle_attr->get<double>("eMin", 0);
        param.estParam.binWidth = wle_attr->get<double>("binWidth", 0.01);
        param.estParam.nBins = wle_attr->get<unsigned int>("nBins", 4040);
    }

    DoS<udrst, feedbackMove, WLEstimator> simulate(rst, rstMove, rnd, param);
    WLEstimator &estm=simulate.getEstimator();
    if (readfile) {
        estm.readHist(readfile);
    }
    simulate.estimate();
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
    return 0;
}
