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
	double value();
	void print_solution(ostream &o) const;
	~rastrigin();
	const double VAR_MAX;
	const double VAR_MIN;
private:
	xmlNode *docroot;
	int dim;
	double *vars;
	unsigned int seed;

};

#endif
