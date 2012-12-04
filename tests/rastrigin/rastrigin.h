#ifndef RASTRIGIN_H
#define RASTRIGIN_H

#include "movable.h"
#include <iostream>
#include <libxml/parser.h>
using namespace std;

class rastrigin
{
public:
	rastrigin(int dimension);
	rastrigin(xmlNode *root);
	int get_dimension() const;
	double get_param(int idx) const;
	void set_param(int idx, double val);
	unsigned int* get_seed() const;
	double value();
	void print_solution(ostream &o) const;
	~rastrigin();
	const double VAR_MAX;
	const double VAR_MIN;
private:
	xmlNode *docroot;
	xmlNode *section;
	int dim;
	double *vars;
	unsigned int seed;

};

#endif
