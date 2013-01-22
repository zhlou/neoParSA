#ifndef RASTRIGIN_H
#define RASTRIGIN_H


#include <iostream>
#include <libxml/tree.h>
using namespace std;
class unirandom;
class rastrigin
{
public:
	rastrigin(int dimension, unirandom &in_rnd);
	rastrigin(xmlNode *root, unirandom &in_rnd);
	int getDimension() const;
	double get_param(int idx) const;
	void set_param(int idx, double val);
	unirandom &rnd;
	unsigned int* get_seed();
	double get_score();
	void generateMove(int idx, double theta_bar);
	void restoreMove(int idx);
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
	double prev_x;
	int prev_idx;
	bool can_rollback;


};

#endif
