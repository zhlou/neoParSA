
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
#include <cmath>

using namespace std;

unirandom::unirandom() : phase(0)
{
	seed = time(NULL);

}

unirandom::unirandom(unsigned int disp) : phase(0)
{
	seed = time(NULL) + disp;
}

unirandom::unirandom(xmlNode* section) : phase(0)
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

// normal(0,1) using rejection sampling a la Knuth TAOCP
double unirandom::randn() 
{
    double x;
    if (0 == phase) {
        do {
            double u1 = random();
            double u2 = random();
            v1 = 2 * u1 - 1;
            v2 = 2 * u2 - 1;
            s = v1 * v1 + v2 * v2;
        } while (s >= 1 || s == 0);
        x = v1 * std::sqrt(-2 * std::log(s) / s);
    } else {
        x = v2 * std::sqrt(-2 * std::log(s) / s);
    }
    phase = 1 - phase;
    return x;
}

double unirandom::laplace(double theta) 
{
    double uniform = 2.0 * random() - 1.0;
    if (uniform >= 0.)
        return -1 * theta * std::log(uniform);
    else
        return theta * std::log(-1*uniform);
}

double unirandom::exponential(double theta) {

    return -1 * theta * std::log(random());
}


double unirandom::lognormal(double mean, double var) 
{
    if (0.0 == var)
        return mean;
    double sigma = std::sqrt(std::log(var/(mean*mean)+1.0));
    double mu = std::log(mean) - sigma * sigma / 2.0;
    return std::exp(randn(mu, sigma));
}


void unirand48::initFromSeed()
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
