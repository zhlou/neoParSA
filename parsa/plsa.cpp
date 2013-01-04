#include "plsa.h"
#include "annealer.h"

plsa::plsa(movable* theproblem, xmlNode* root, MPI_Comm thecomm, int in_nnodes,
        int in_rank) :
        lam(theproblem, root), comm(thecomm), nnodes(in_nnodes), rank(in_rank),
        StatsComm(thecomm, in_nnodes, in_rank)
{
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
    StatsComm.CommSegment(mean, vari, fit_mean, fit_sd);


}
