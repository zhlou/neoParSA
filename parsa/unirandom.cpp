
/*
 * unirandom.cpp
 *
 *  Created on: Dec 5, 2012
 *      Author: zhlou
 */

#include "unirandom.h"
#include <cstdlib>
#include <ctime>

using namespace std;

unirandom::unirandom()
{
	seed = time(NULL);

}

unirandom::unirandom(unsigned int in_seed)
{
	seed = in_seed;
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
