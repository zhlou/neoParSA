/*
 * normSchwefel.h
 *
 *  Created on: Oct 12, 2015
 *      Author: zhlou
 */

#ifndef TESTS_SCHWEFEL_NORMSCHWEFEL_H_
#define TESTS_SCHWEFEL_NORMSCHWEFEL_H_

#include <iostream>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
class unirandom;

class normSchwefel {
public:
	normSchwefel(int dimension, unirandom &in_rnd);
	//normSchwefel(xmlNode *root, unirandom &in_rnd, const char *secName="Schwefel");
    normSchwefel(ptree &root, unirandom &in_rnd, std::string secName="Schwefel");
	~normSchwefel();
	int getDimension() const {return 1;}
	double getParam(int idx) const;
	void setParam(int idx, double val);
	unirandom &rnd;
	double get_score();
	void generateMove(int, double theta);
	void restoreMove(int);
	void printSolution(std::ostream &o) const;
	//void writeSection(xmlNode *docroot, xmlChar *secname);
    void writeSection(ptree &root, std::string secname);

	int getStateSize() const {return sizeof(double)*dim;}
	void serialize(void *buf) const;
	void deserialize(void const *buf);
	double scramble();

	static const double VAR_MAX;
	static const double VAR_MIN;

private:
	//xmlNode *docroot;
	//xmlNode *section;
	int dim;
	int idx;
	double *vars;
	double prev_x;
	int prev_idx;
	bool outOfBounds;
	bool can_rollback;
};

#endif /* TESTS_SCHWEFEL_NORMSCHWEFEL_H_ */
