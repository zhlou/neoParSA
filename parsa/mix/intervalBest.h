/* intervalBest.h */
#ifndef INTERVALBEST_H
#define INTERVALBEST_H

#include "dynDebug.h"
#include "MPIState.h"
#include "mixState.h"
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

template <class Problem>
class intervalBest {
public:
	class Param {
	public:
		int interval;
        Param(const ptree &root);
	};
	intervalBest(Problem &in_problem, const MPIState &mpiState,
			unirandom &in_rand, const Param &param);
	~intervalBest();
	mixState Mix(aState &state);
	static const char * name;
	void setDebug(debugStatus st, const char* outname=NULL){
		debugOut.setDebug(ignore, NULL);
	}
private:
	Problem &problem;
	const int interval;
	const MPIState &mpi;
	unirandom &rnd;
	dynDebug debugOut;
	int count;
	int buf_size;
	void *state_buf;
};


#include "intervalBest.hpp"

#endif
