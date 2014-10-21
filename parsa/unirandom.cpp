
/*
 * unirandom.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: zhlou
 */

#include "unirandom.h"
#include <cstdlib>
#include <ctime>
#include <climits>

using namespace std;

unirandom::unirandom()
{
	seed = time(NULL);

}

unirandom::unirandom(unsigned int disp)
{
	seed = time(NULL) + disp;
}

unirandom::unirandom(xmlNode* section)
{
	xmlChar *prop = NULL;
	if ((prop = xmlGetProp(section, (xmlChar *)"seed")) != NULL) {
		seed = atoi((char * )prop);
		xmlFree(prop);
		prop = NULL;
	} else
		seed = time(NULL);
}

double unirandom::random()
{
	return (double)rand_r(&seed)/RAND_MAX;
}
inline void unirand48::initFromSeed()
{
    xsubi[0] = 0x330E; // follow what srand48 do
    xsubi[1] = static_cast<unsigned short>(seed);
    xsubi[2] = static_cast<unsigned short>(seed >> (CHAR_BIT * sizeof(unsigned short)));
}
unirand48::unirand48()
{
    initFromSeed();
}

unirand48::unirand48(unsigned int disp) : unirandom(disp)
{
    initFromSeed();
}

unirand48::unirand48(xmlNode *section) : unirandom(section)
{
    initFromSeed();
}
