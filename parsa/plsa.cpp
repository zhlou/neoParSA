#include "plsa.h"
#include "annealer.h"

plsa::plsa(movable* theproblem, xmlNode* root, MPI_Comm thecomm) :
        lam(theproblem, root), comm(thecomm)
{
    MPI_Comm_size(comm, &nsize);
    MPI_Comm_rank(comm, &rank);
}

plsa::~plsa()
{

}

/*
double plsa::loop()
{
    long unsigned step_cnt = 0;
    int accept = 0, i;
    double delta, vari;
    while (!frozen()) {
    	accept = 0;
    	vari = 0.;
    	for (i = 0; i < proc_tau; i++) {
    		if ((delta = move()) != 0.)
    			accept ++;

    	}

    }

    return problem->get_score();
}
 */

void plsa::updateSegment()
{
    lam::updateSegment();

}
