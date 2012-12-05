/*
 * unirandom.h
 *
 *  Created on: Dec 5, 2012
 *      Author: zhlou
 */

#ifndef UNIRANDOM_H_
#define UNIRANDOM_H_
#include <libxml/tree.h>

/*
 * This class encapsulates the generation of a random value
 * from uniform [0,1].
 */
class unirandom
{
public:
	unirandom();
	unirandom(unsigned int in_seed);
	unirandom(xmlNode *section);
	double random();
private:
	unsigned int seed;
};


#endif /* UNIRANDOM_H_ */
