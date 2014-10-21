#ifndef UDRST_H
#define UDRST_H


#include <iostream>
#include <libxml/tree.h>
using namespace std;
class unirandom;
// Rastrigin function with uniform dimension move generation
class udrst
{
public:
	udrst(int dimension, unirandom &in_rnd);
	udrst(xmlNode *root, unirandom &in_rnd);
	int getDimension() const {return 1;};
	double get_param(int idx) const;
	void set_param(int idx, double val);
	unirandom &rnd;
	//unsigned int* get_seed();
	double get_score();
	void generateMove(int, double theta);
	void restoreMove(int);
	void print_solution(ostream &o) const;
	void write_section(xmlChar *secname);
	~udrst();
	static const double VAR_MAX;
	static const double VAR_MIN;

	int getStateSize() const {return sizeof(double)*dim;};
	void serialize(void *buf) const;
	void deserialize(void const *buf);
	double scramble();
private:
	xmlNode *docroot;
	xmlNode *section;
    int dim;
    int idx;
	double *vars;
	double prev_x;
	int prev_idx;
	bool outOfBounds;
	bool can_rollback;


};

#endif
