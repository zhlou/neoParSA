#include "plsa.h"
#include <cmath>
using namespace std;

const char *plsa::name = "parallel Lam";

plsa::Param::Param(xmlNode *root, debugStatus st, const char *outname) :
        lamParam(root, st, outname)
{

}

plsa::plsa(Param param, const MPIState &mpiState) :
        lam(param.lamParam), mpi(mpiState)
{
    local_stat_buf = new StatData[mpi.nnodes];
    //l_stat.s = -1;
    lambda *= (double)mpi.nnodes;
    MPI_Win_create(&l_stat, sizeof(StatData), sizeof(StatData), MPI_INFO_NULL,
            mpi.comm, &stat_win);
}

plsa::~plsa()
{
    MPI_Win_free(&stat_win);
    delete[] local_stat_buf;
}

/*

bool plsa::global_frozen()
{
    int local_freeze = (freeze_cnt >= cnt_crit);
    int global_flag;
    MPI_Allreduce(&local_freeze, &global_flag, 1, MPI_INT, MPI_SUM, mpi.comm);
    return (bool)global_flag;
}

 */

/* This method takes a boolean to decide the content of l_stat.var
 * True (default) for using standard deviation or square root of the value
 * of vari. False for using vari itself. */
void plsa::PackNCommStats(bool UseSD)
{
    //pack statistics
    //l_stat.s = s;
    l_stat.mean = mean;
    if (UseSD)
        l_stat.var = sqrt(vari);
    else
        l_stat.var = vari;
    l_stat.success = success;
    //l_stat.moves = proc_tau;
    //l_stat.energy = energy;
    local_stat_buf[mpi.rank] = l_stat;
    MPI_Win_post(mpi.group, 0, stat_win);
    MPI_Win_start(mpi.group, 0, stat_win);
    for (int j = 0; j < mpi.nnodes; ++j) {
        if (j != mpi.rank) {
            MPI_Get(local_stat_buf + j, sizeof(StatData), MPI_BYTE, j, 0,
                    sizeof(StatData), MPI_BYTE, stat_win);
        }
    }
    MPI_Win_complete(stat_win);
    MPI_Win_wait(stat_win);
}

void plsa::updateEstimators(double s)
{
    PackNCommStats();
    // local_stat_buf is not ready for reading until all the communications are
    // finished. So two separate loops are necessary.
    double w_m = 1.0 / (double)mpi.nnodes;
    fit_mean->decay();
    fit_sd->decay();
    for (int i = 0; i < mpi.nnodes; ++i) {
        fit_mean->partialUpdate(w_m, 1.0 / local_stat_buf[i].mean, s);
        fit_sd->partialUpdate(w_m, 1.0 / local_stat_buf[i].var, s);
        success += local_stat_buf[i].success;
    }
    fit_mean->finishUpdate();
    fit_sd->finishUpdate();
    acc_ratio = (double) success / (double)(mpi.nnodes * proc_tau);
}

void plsa::collectInitStats(unsigned long init_loop)
{
    PackNCommStats(false);
    mean = 0;
    vari = 0;
    success = 0;
    for (int i = 0; i < mpi.nnodes; ++i) {
        mean += local_stat_buf[i].mean;
        vari += local_stat_buf[i].var;
        success += local_stat_buf[i].success;
    }
    double total_moves = (double)mpi.nnodes * init_loop;
    acc_ratio = (double) success / total_moves;
    mean /= total_moves;
    vari = vari / total_moves - mean * mean;
}

void plsa::collectInitStats(double initMean, double initVar, double initAccRatio) 
{
    double var[3] = {initMean, initVar, initAccRatio};
    MPI_Allreduce(MPI_IN_PLACE, var, 3, MPI_DOUBLE, MPI_SUM, mpi.comm);
    mean = var[0] / mpi.nnodes;
    vari = var[1] / mpi.nnodes;
    acc_ratio = var[2] / mpi.nnodes;   
    
}
