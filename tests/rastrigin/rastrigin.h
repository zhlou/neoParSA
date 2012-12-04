#ifndef RASTRIGIN_H
#define RASTRIGIN_H

#include "movable.h"
#include <iostream>
#include <libxml/tree.h>
using namespace std;

class rastrigin
{
public:
	rastrigin(int dimension);
	rastrigin(xmlNode *root);
	int get_dimension() const;
	double get_param(int idx) const;
	void set_param(int idx, double val);
	unsigned int* get_seed();
	double value();
	void print_solution(ostream &o) const;
	void write_section(xmlChar *secname);
	~rastrigin();
	static const double VAR_MAX;
	static const double VAR_MIN;
	int dim;
private:
	xmlNode *docroot;
	xmlNode *section;
	double *vars;
	unsigned int seed;

};

#endif
