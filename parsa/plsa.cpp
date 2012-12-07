#include "plsa.h"
#include "annealer.h"

plsa::plsa(movable* theproblem, xmlNode* root, MPI_Comm thecomm) :
        annealer(theproblem, root), comm(thecomm)
{
    MPI_Comm_size(comm, &nsize);
    MPI_Comm_rank(comm, &rank);
    proc_tau = 100; // for now

}

plsa::~plsa()
{

}

double plsa::loop()
{
    long unsigned step_cnt = 0;
    while (!frozen()) {

    }

    return problem->get_score();
}
