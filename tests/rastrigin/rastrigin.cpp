#include <cmath>
#include <cstdlib>
#include <ctime>
#include <libxml/tree.h>
#include <string>
#include <sstream>
#include "rastrigin.h"

using namespace std;
const double rastrigin::VAR_MAX = 5.12;
const double rastrigin::VAR_MIN = -5.12;

rastrigin::rastrigin(int dimension) :
		dim(dimension)
{
	int i;
	seed = time(NULL);
	vars = new double[dim];
	for (i = 0; i < dim; i++) {
		vars[i] = VAR_MAX * ((rand_r(&seed) * 2.0 / RAND_MAX) - 1.0);
	}
	docroot = NULL;
	section = NULL;
}

rastrigin::rastrigin(xmlNode *root)
{
	docroot = root;
	section = root->children;
	xmlChar *prop = NULL;
	while (section != NULL) {
		if (!xmlStrcmp(section->name, (xmlChar *) "rastrigin"))
			break;
		section = section->next;
	}
	if (section == NULL) {
		throw 1;
	}
	if ((prop = xmlGetProp(section, (xmlChar *) "dim")) != NULL) {
		dim = atoi((char *) prop);
		xmlFree(prop);
		prop = NULL;
	} else {
		throw 2;
	}
	vars = new double[dim];
	if ((prop = xmlGetProp(section, (xmlChar *) "seed")) != NULL) {
		seed = atoi((char *) prop);
		xmlFree(prop);
		prop = NULL;
	} else {
		seed = time(NULL);
	}
	char *namebuf = new char[255];
	for (int i = 0; i < dim; i++) {
		sprintf(namebuf, "x%d", i+1);
		if ((prop = xmlGetProp(section, (xmlChar *)namebuf))
				!= NULL) {
			vars[i] = strtod((char *) prop, NULL);
			xmlFree(prop);
			prop = NULL;
		} else {
			vars[i] = VAR_MAX * ((rand_r(&seed) * 2.0 / RAND_MAX) - 1.0);
		}

	}

}

void rastrigin::write_section(xmlChar *secname)
{
	xmlNode *node = xmlNewChild(docroot, NULL, secname, NULL);
	char *namebuf = new char[255];
	char *valbuf = new char[255];
	sprintf(namebuf, "%d", dim);
	xmlNewProp(node, (xmlChar *) "dim", (xmlChar *) namebuf);
	for (int i = 0; i < dim; i++) {
		sprintf(namebuf, "x%d", i + 1);
		sprintf(valbuf, "%f", vars[i]);
		xmlNewProp(node, (xmlChar *) namebuf, (xmlChar *) valbuf);
	}
	delete[] valbuf;
	delete[] namebuf;

}

void rastrigin::print_solution(ostream& o) const
{
	o << "{" << endl;
	for (int i = 0; i < dim; i++) {
		o << vars[i] << endl;
	}
	o << "}" << endl;
}

rastrigin::~rastrigin()
{
	if (vars != NULL) {
		delete[] vars;
	}
	vars = NULL;
}
double rastrigin::value()
{
	double tot = 0;
	for (int i = 0; i < dim; i++) {
		tot += vars[i] * vars[i] - 10.0 * cos(2 * M_PI * vars[i]);
	}
	return (10 * dim + tot);
}

int rastrigin::get_dimension() const
{
	return dim;
}

double rastrigin::get_param(int idx) const
{
		return vars[idx];
}

void rastrigin::set_param(int idx, double val)
{
	vars[idx] = val;
}

unsigned int* rastrigin::get_seed()
{
	return &seed;
}
