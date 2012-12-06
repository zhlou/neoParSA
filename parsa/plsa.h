/*
 * plsa.h
 *
 *  Created on: Dec 6, 2012
 *      Author: zhlou
 */

#ifndef PLSA_H_
#define PLSA_H_
#include <mpi.h>
#include <libxml/tree.h>
#include "annealer.h"


class plsa : public annealer
{
public:
	plsa(movable *theproblem, xmlNode *root, MPI_Comm thecomm);
	~plsa();
protected:
	void cool_s();
	bool frozen();
	MPI_Comm comm;
private:

};


#endif /* PLSA_H_ */
