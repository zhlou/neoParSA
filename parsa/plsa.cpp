#include "plsa.h"


plsa::plsa(movable* theproblem, xmlNode* root, MPI_Comm thecomm, int in_nnodes,
        int in_rank) :
        lam(theproblem, root), comm(thecomm), nnodes(in_nnodes), rank(in_rank)
{
    local_stat_buf = new StatData[nnodes];
    l_stat.s = -1;
    MPI_Win_create(&l_stat, sizeof(StatData), sizeof(StatData), MPI_INFO_NULL,
            comm, &stat_win);
    MPI_Comm_group(comm, &group);
}

plsa::~plsa()
{
    MPI_Win_free(&stat_win);
    delete[] local_stat_buf;
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
    collectStats();
    CommSegment();
    updateLam();
}

void plsa::CommSegment()
{
    l_stat.s = s;
    l_stat.mean = mean;
    l_stat.sd = sqrt(vari);
    l_stat.success = success;
    //l_stat.moves = proc_tau;
    l_stat.energy = energy;
    local_stat_buf[rank] = l_stat;
    MPI_Win_post(group, 0, stat_win);
    MPI_Win_start(group, 0, stat_win);
    for (int j = 0; j < nnodes; ++j) {
        if (j != rank) {
            MPI_Get(local_stat_buf + j, sizeof(StatData), MPI_BYTE, j, 0,
                    sizeof(StatData), MPI_BYTE, stat_win);
        }

    }
    MPI_Win_complete(stat_win);
    MPI_Win_wait(stat_win);
    double w_m = 1.0 / (double)nnodes;
    for (int i = 0; i < nnodes; ++i) {
        fit_mean(w_m, s, local_stat_buf[i].mean);
        fit_sd(w_m, s, local_stat_buf[i].sd);
        success += local_stat_buf[i].success;
    }
    acc_ratio = (double) success / (double)(nnodes * proc_tau);
}
