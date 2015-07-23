/*
 * rstWLS.cpp
 */
#include <iostream>
#include<libxml/parser.h>
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
    char *docname = argv[1];
    char *outname = NULL;
    if (argc == 3)
        outname = argv[2];
    xmlDocPtr xmldoc = xmlParseFile(docname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);
    
    rastrigin rst(xmlroot, rnd);
    feedbackMove<rastrigin> rstMove(rst, rnd, xmlroot);
    DoS<rastrigin, feedbackMove, WLEstimator>::Param param;
    param.initWeight = 1e-2;
    param.nSteps = 100000;
    param.estParam.eMin=0;
    param.estParam.binWidth=0.01;
    param.estParam.nBins=4040;
    xmlNodePtr DoSParamNode=getSectionByName(xmlroot,"DoS");
    if (DoSParamNode != NULL) {
        param.initWeight = getPropDouble(DoSParamNode, "weight");
        param.nSteps = getPropLong(DoSParamNode,"nsteps");
    }
    xmlNodePtr WLENode=getSectionByName(xmlroot, "WLEstimator");
    if (WLENode != NULL) {
        param.estParam.eMin = getPropDouble(WLENode,"eMin");
        param.estParam.binWidth = getPropDouble(WLENode,"binWidth");
        param.estParam.nBins = getPropInt(WLENode, "nBins");
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