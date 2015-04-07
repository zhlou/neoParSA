/* 
 * File:   coupledAnnealer.hpp
 * Author: zhlou
 *
 * Created on January 28, 2015, 5:32 PM
 */

#ifndef COUPLEDANNEALER_HPP
#define	COUPLEDANNEALER_HPP



template <class Problem, class Schedule, class FrozenCnd, 
        template<class> class CoupledMove>
coupledAnnealer<Problem, Schedule, FrozenCnd, CoupledMove>::coupledAnnealer(
        Problem& problem, 
        unirandom& in_rand, 
        typename Schedule::Param scheParam, 
        typename FrozenCnd::Param frozenParam, 
        const typename CoupledMove<Problem>::Param& moveParam, 
        xmlNode* root, 
        const MPIState& mpiState) :
        annealer<Problem, Schedule, FrozenCnd, CoupledMove>::annealer(problem, in_rand,root),
        mpi(mpiState)
{
    this->cooling = new Schedule(scheParam, mpiState);
    this->move = new CoupledMove<Problem>(problem, mpiState, in_rand, moveParam);
    this->frozen = new FrozenCnd(frozenParam, mpiState);
    this->state.energy = this->move->get_score();
}

template<class Problem, class Schedule, class FrozenCnd, 
        template<class> class CoupledMove>
coupledAnnealer<Problem, Schedule, FrozenCnd, CoupledMove>::~coupledAnnealer()
{
    
}

template<class Problem, class Schedule, class FrozenCnd, 
        template<class> class CoupledMove>
int coupledAnnealer<Problem, Schedule, FrozenCnd, CoupledMove>::getWinner()
{
    struct
    {
        double energy;
        int rank;
    } doubleint;
    doubleint.energy = (this->state).energy;
    doubleint.rank = mpi.rank;
    MPI_Allreduce(MPI_IN_PLACE, &doubleint, 1, MPI_DOUBLE_INT, MPI_MINLOC,
            mpi.comm);
    return doubleint.rank;
}

template<class Problem, class Schedule, class FrozenCnd, 
        template<class> class CoupledMove>
void coupledAnnealer<Problem, Schedule, FrozenCnd, CoupledMove>::updateStats(aState& state)
{
    annealer<Problem, Schedule, FrozenCnd, CoupledMove>::updateStats(state);
    this->move->Mix(state);
}

#endif	/* COUPLEDANNEALER_HPP */

