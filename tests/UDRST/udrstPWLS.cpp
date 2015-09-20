/*
 * rstWLS.cpp
 */
#include <iostream>
#include <exception>
#include <libxml/parser.h>
#include <unistd.h>
#include "DoS/DoS.h"
#include "DoS/PWLE.h"
#include "move/feedbackMove.h"
#include "udrst.h"
#include "unirandom.h"




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

    char *docname = argv[optind];
    xmlDocPtr xmldoc = xmlParseFile(docname);
    xmlNodePtr xmlroot = xmlDocGetRootElement(xmldoc);
    unirandom rnd(mpi.rank);
    udrst rst(xmlroot, rnd);
    feedbackMove<udrst> rstMove(rst, rnd, xmlroot);
    DoS<udrst, feedbackMove, PWLE>::Param param;
    param.initWeight = 1e-2;
    param.nSteps = 100000;
    param.estParam.eMin=0;
    param.estParam.binWidth=0.01;
    param.estParam.nBins=4040;
    param.estParam.mpi = &mpi;
    param.estParam.syncFreq = 1000;
    param.estParam.saveName = NULL;
    param.estParam.saveFreq = 0;
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
        param.estParam.syncFreq = getPropInt(WLENode, "syncFreq");
        try {
            param.estParam.saveFreq = getPropInt(WLENode, "saveFreq");
        } catch (const std::exception &e) {
            // ignored
        }
        param.estParam.saveName = (char *)xmlGetProp(WLENode, BAD_CAST"saveName");
    }
    DoS<udrst, feedbackMove, PWLE> simulate(rst, rstMove, rnd, param);
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

    if (param.estParam.saveName) {
        xmlFree(param.estParam.saveName);
    }
    xmlFreeDoc(xmldoc);
    xmlCleanupParser();
    MPI_Finalize();
    return 0;
}
