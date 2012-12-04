#include <cmath>
#include <cstdlib>
#include <ctime>
#include <libxml/tree.h>
#include <string>
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
	if ((prop = xmlGetProp(section, (xmlChar *)"dim")) != NULL) {
		dim = atoi((char *)prop);
		xmlFree(prop);
		prop = NULL;
	} else {
		throw 2;
	}
	vars = new double[dim];
	if ((prop = xmlGetProp(section, (xmlChar *)"seed")) != NULL) {
		seed = atoi((char *)prop);
		xmlFree(prop);
		prop = NULL;
	} else {
		seed = time(NULL);
	}
	string str("x");
	for (int i = 0; i < dim; i++) {
		if ((prop = xmlGetProp(section, (xmlChar *)(str + (i+1)).c_str()))
				!= NULL) {
			vars[i] = strtod((char *)prop, NULL);
			xmlFree(prop);
			prop = NULL;
		} else {
			vars[i] = VAR_MAX * ((rand_r(&seed) * 2.0 / RAND_MAX) - 1.0);
		}

	}

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

int rastrigin::get_dimension() const {
	return dim;
}
