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
#include "annealer.h"
#include "parallelFBMove.h"
#include "unirandom.h"
#include "plsa.h"
#include "debugOut.h"
using namespace std;
int main(int argc, char **argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
    unirandom rnd(rank);
    rastrigin rst(docroot, rnd);
    parallelFBMove<rastrigin, debugIGNORE> rst_problem(rst, docroot, MPI_COMM_WORLD,
            size, rank);
    plsa *pschedule = new plsa(docroot, MPI_COMM_WORLD, size, rank);

    //feedbackMove<rastrigin> rst_problem(rst, docroot);
    //lam schedule(docroot);
    annealer<plsa, parallelFBMove<rastrigin, debugIGNORE>, unirandom >
            rst_anneal(*pschedule, rst_problem, rnd, docroot);
    cout << "The initial state is " << endl;
    rst.print_solution(cout);
    cout << "The fininal energy is " << rst_anneal.loop() << endl;
    cout << "The solution is " << endl;
    rst.print_solution(cout);
    delete pschedule;
    MPI_Finalize();

    return 0;
}

