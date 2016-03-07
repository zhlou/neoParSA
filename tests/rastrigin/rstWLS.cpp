/*
 * rstWLS.cpp
 */
#include <iostream>
#include <string>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#include <boost/optional.hpp>

#include "DoS/DoS.h"
#include "DoS/WLEstimator.h"
#include "move/feedbackMove.h"
#include "rastrigin.h"
#include "unirandom.h"




int main(int argc, char **argv)
{
    if (argc <= 1) {
        std::cerr << "Missing input files" << std::endl;
        std::cerr << "Usage: rstWLS input [output]" << std::endl;
        return 1;
    }
    unirandom rnd;
    std::string docname(argv[1]);
    char *outname = NULL;
    if (argc == 3)
        outname = argv[2];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &xmlroot = pt.begin()->second;

    rastrigin rst(xmlroot, rnd);
    feedbackMove<rastrigin> rstMove(rst, rnd, xmlroot);
    DoS<rastrigin, feedbackMove, WLEstimator>::Param param;
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

    DoS<rastrigin, feedbackMove, WLEstimator> simulate(rst, rstMove, rnd, param);
    simulate.estimate();
    WLEstimator &estm=simulate.getEstimator();
    if (outname == NULL)
        estm.printHist(std::cout);
    else {
        ofstream myout(outname);
        estm.printHist(myout);
        myout.close();
    }
    return 0;
}
