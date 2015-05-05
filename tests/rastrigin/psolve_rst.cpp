/*
 * psolve_rst.cpp
 *
 *  Created on: Jan 28, 2013
 *      Author: zhlou
 */



#include <iostream>
#include <libxml/parser.h>
#include <mpi.h>
#include "rastrigin.h"
#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "unirandom.h"
#include "plsa.h"
#include "dynDebug.h"
#include "adaptMix.h"
#include "criCountP.h"
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
        cerr<< "Missing input files" << endl;
        return 1;
    }
    char *docname = argv[1];
    xmlDoc *doc = xmlParseFile(docname);
    xmlNode *docroot = xmlDocGetRootElement(doc);
    if (doc == NULL) {
        cerr << "Input incorrect" << endl;
        return 2;
    }
    unirandom rnd(mpi.rank);
    rastrigin rst(docroot, rnd);
    //parallelFBMove<rastrigin, debugIGNORE, adaptMix> *rst_problem =
    //        new parallelFBMove<rastrigin, debugIGNORE, adaptMix>(rst, rnd,
    //                                                             docroot, mpi);
    //plsa *pschedule = new plsa(docroot, mpi);
    plsa::Param scheduleParam(docroot);
    criCountP::Param criCntParam(docroot);
    adaptMix<rastrigin>::Param mixParam(docroot);

    //feedbackMove<rastrigin> rst_problem(rst, docroot);
    //lam schedule(docroot);
    //annealer<plsa, parallelFBMove<rastrigin, debugIGNORE, adaptMix>, unirandom >
    //        rst_anneal(*pschedule, *rst_problem, rnd, docroot);
    pannealer<rastrigin, plsa, criCountP, parallelFBMove, adaptMix> *rst_anneal =
            new pannealer<rastrigin, plsa, criCountP, parallelFBMove, adaptMix>
            (rst, rnd, scheduleParam, criCntParam, mixParam, docroot, mpi);
    cout << "The initial state is " << endl;
    rst.print_solution(cout);
    cout << "The fininal energy is " << rst_anneal->loop() << endl;
    cout << "The solution is " << endl;
    rst.print_solution(cout);
    if (rst_anneal->getWinner() == mpi.rank) {
        cout << "The winning energy is " << rst.get_score() << endl;
        cout << "The winning solution is " << endl;
        rst.print_solution(cout);
    }
    //delete pschedule;
    delete rst_anneal;
    MPI_Finalize();

    return 0;
}

