/* 
 * File:   FBMoveIntervalMix.hpp
 * Author: zhlou
 *
 * Created on January 30, 2015, 2:02 PM
 */

#ifndef FBMOVEINTERVALMIX_HPP
#define	FBMOVEINTERVALMIX_HPP

#include <cmath>
#include <limits>


using namespace std;

template<class Problem>
const char * FBMoveIntervalMix<Problem>::name = "FeedbackMove+IntervalMix";

template<class Problem>
FBMoveIntervalMix<Problem>::Param::Param() : mix_interval(100),
        move_interval(100), move_gain(0.03), target(0.44), initTheta(1.0),
        thetaMin(0.), thetaMax(numeric_limits<double>::max()), mix_target(0.5)
{

}

template<class Problem>
FBMoveIntervalMix<Problem>::Param::Param(xmlNode *root) : mix_interval(100),
        move_interval(100), move_gain(0.03), target(0.44), initTheta(1.0),
        thetaMin(0.), thetaMax(numeric_limits<double>::max()), mix_target(0.5)
{
    if (!root) {
        xmlNode *section = getSectionByName(root, "moveMix");
        if (!section) {
            try {
                mix_interval = getPropInt(section, "mix_interval");
            } catch (const std::exception &e) {

            }
            try {
                move_interval = getPropInt(section, "move_interval");
            } catch (const std::exception &e) {

            }
            try {
                move_gain = getPropDouble(section, "move_gain");
            } catch (const std::exception &e) {

            }
            try {
                target = getPropDouble(section, "move_target");
            } catch (const std::exception &e) {

            }
            try {
                initTheta = getPropDouble(section, "initTheta");
            } catch (const std::exception &e) {

            }
            try {
                thetaMin = getPropDouble(section, "thetaMin");
            } catch (const std::exception &e) {

            }
            try {
                thetaMax = getPropDouble(section, "thetaMax");
            } catch (const std::exception &e) {

            }
            try {
                mix_target = getPropDouble(section, "mix_target");
            } catch (const std::exception &e) {

            }

        }
    }

}

template<class Problem>
FBMoveIntervalMix<Problem>::FBMoveIntervalMix(Problem& in_problem,
        const MPIState& mpiState,
        unirandom& in_rnd,
        const Param& param) :
        problem(in_problem),
        rnd(in_rnd),
        mpi(mpiState),
        nparams(problem.getDimension()),
        moveCore(nparams, param.move_gain, param.target, param.initTheta, param.thetaMin,
        param.thetaMax, mpi, rnd, debugOut),
        mix_interval(param.mix_interval),
        move_interval(param.move_interval),
        index(-1),
        sweep(0),
        tau_count(0),
        mix(in_problem, mpiState, in_rnd),
        mix_target(param.mix_target)
{
    energy = problem.get_score();
    prev_energy = energy;
    adoptArray = new int[mpi.nnodes];

}

template<class Problem>
FBMoveIntervalMix<Problem>::~FBMoveIntervalMix()
{
    delete[] adoptArray;
}

template<class Problem>
double FBMoveIntervalMix<Problem>::propose()
{
    index++;
    index %= nparams;
    if (index == 0) {
        ++sweep;
        if (sweep % move_interval == 0) {
            move_control();
        }
    }

    prev_energy = energy;

    problem.generateMove(index, moveCore.genMove(index));
    energy = problem.get_score();
    return (energy - prev_energy);
}

template<class Problem>
void FBMoveIntervalMix<Problem>::accept()
{
    moveCore.accept(index);
}

template<class Problem>
void FBMoveIntervalMix<Problem>::reject()
{
    moveCore.reject(index);
    problem.restoreMove(index);
    energy = prev_energy;
}

template<class Problem>
void FBMoveIntervalMix<Problem>::move_control()
{
    debugOut << sweep;
    moveCore.moveControl();

    if (!debugOut.isIgnore()) {
        for (int i = 0; i < nparams; ++i) {
            debugOut << "\t" << moveCore.getThetaBar(i);
        }
    }

    debugOut << endl;

}

template<class Problem>
mixState FBMoveIntervalMix<Problem>::Mix(aState& state)
{
    ++tau_count;
    if ((tau_count % mix_interval) != 0)
        return mixState();
    int i, p, nadopt = 0;
    tau_count = 0;
    mix.calProbTab(state);
    p = mix.getPartner();
    state.energy = mix.adoptState(i);
    for (i = 0; i < mpi.nnodes; ++i) {
        adoptArray[i] = 0;
    }
    adoptArray[p] = 1;
    MPI_Allreduce(MPI_IN_PLACE, adoptArray, mpi.nnodes, MPI_INT, MPI_LOR,
            mpi.comm);
    for (i = 0; i < mpi.nnodes; ++i)
        nadopt += adoptArray[i];
    double adoptRate = (double) (nadopt) / mpi.nnodes;
    if (adoptRate < mix_target) {
        moveCore.setVarTheta(1.0 / adoptRate - 1 / mix_target);
    } else {
        moveCore.setVarTheta(0.);
    }


    return mixState(p);

}

#endif	/* FBMOVEINTERVALMIX_HPP */

