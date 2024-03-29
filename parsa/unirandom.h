/*
 * unirandom.h
 *
 *  Created on: Dec 5, 2012
 *      Author: zhlou
 */

#ifndef UNIRANDOM_H_
#define UNIRANDOM_H_
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;

#include <cstdlib>
using namespace std;
/*
 * This class encapsulates the generation of a random value
 * from uniform [0,1).
 */
class unirandom
{
public:
    virtual ~unirandom() {}
    unirandom();
    unirandom(unsigned int disp);
    unirandom(const ptree &root);
    virtual double random();
    double randn();
    double randn(double mu, double sigma) {return sigma*randn()+mu;}
    double laplace(double theta);
    double lognormal(double mean, double var);
    double exponential(double theta);
    double cauchy(double gamma); // symmetric cauchy
    virtual void setSeed(unsigned int newSeed)
    {
        seed = newSeed;
    }
protected:
    unsigned int seed;
private:
    double v1, v2, s;
    int phase;
};

class unirand48: public unirandom
{
private:
    unsigned short xsubi[3];
    void initFromSeed();
public:
    unirand48();
    unirand48(unsigned int disp);
    unirand48(ptree &ptree);
    void setSeed(unsigned int newSeed)
    {
        seed = newSeed;
        initFromSeed();
    }
    virtual double random(){return erand48(xsubi);}

};

#endif /* UNIRANDOM_H_ */
