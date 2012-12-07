#include "plsa.h"
#include "annealer.h"

plsa::plsa(movable* theproblem, xmlNode* root, MPI_Comm thecomm) : 
annealer(theproblem, root), comm(thecomm)
{
    MPI_Comm_size(comm, &nsize);
    MPI_Comm_rank(comm, &rank);

}

plsa::~plsa()
{
    
}

void plsa::cool_s()
{
    
}