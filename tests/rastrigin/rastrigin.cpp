#include <cmath>
#include <cstdlib>
#include <ctime>
#include <libxml/tree.h>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>
#include "rastrigin.h"
#include "unirandom.h"
#include "xmlUtils.h"
#include <string.h>
// #include <omp.h>

using namespace std;
const double rastrigin::VAR_MAX = 5.12;
const double rastrigin::VAR_MIN = -5.12;

rastrigin::rastrigin(int dimension, unirandom &in_rnd) :
		dim(dimension), rnd(in_rnd)
{
	int i;
	vars = new double[dim];
	for (i = 0; i < dim; i++) {
		vars[i] = VAR_MAX * (rnd.random() *2.0 -1.0);
	}
	docroot = NULL;
	section = NULL;
	prev_x = 0; // to make the compiler happy
	prev_idx = -1;
	can_rollback = false;
	outOfBounds = false;

}

rastrigin::rastrigin(xmlNode *root, unirandom &in_rnd):
		rnd(in_rnd)
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
	char *namebuf = NULL;
	for (int i = 0; i < dim; i++) {
		asprintf(&namebuf, "x%d", i+1);
		if ((prop = xmlGetProp(section, (xmlChar *)namebuf))
				!= NULL) {
			vars[i] = strtod((char *) prop, NULL);
			xmlFree(prop);
			prop = NULL;
		} else {
			vars[i] = VAR_MAX * (2.0 * rnd.random() - 1.0);
			// write the value back into xml file
			char *valbuf = NULL;
			asprintf(&valbuf, "%g", vars[i]);
			xmlNewProp(section, (const xmlChar *)namebuf,
			        (const xmlChar*)valbuf);
			free(valbuf);
		}
		free(namebuf);
		namebuf = NULL;

	}


    prev_x = 0; // to make the compiler happy
    prev_idx = -1;
    can_rollback = false;
    outOfBounds = false;
}

void rastrigin::write_section(xmlChar *secname)
{
    xmlNode *node;
    node = getSectionByName(docroot, (const char *)secname);
    if (node != NULL) {
        xmlUnlinkNode(node);
        xmlFreeNode(node);
    }
	node = xmlNewChild(docroot, NULL, secname, NULL);
	char *namebuf = new char[255];
	char *valbuf = new char[255];
	sprintf(namebuf, "%d", dim);
	xmlNewProp(node, (xmlChar *) "dim", (xmlChar *) namebuf);
	for (int i = 0; i < dim; i++) {
		sprintf(namebuf, "x%d", i + 1);
		sprintf(valbuf, "%g", vars[i]);
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

void rastrigin::generateMove(int idx, double theta)
{
    double x = prev_x = vars[idx];
    x += theta;
    if (x > VAR_MAX || x < VAR_MIN) {
        x = prev_x;
        outOfBounds = true;
    }
//    x = fmod(x,VAR_MAX);

//    while (x > rastrigin::VAR_MAX || x < rastrigin::VAR_MIN) {
//        if (x > rastrigin::VAR_MAX)
//            x = 2 * rastrigin::VAR_MAX - x;
//        if (x < rastrigin::VAR_MIN)
//            x = 2 * rastrigin::VAR_MIN - x;
//    }
    vars[idx] = x;
    can_rollback = true;
    prev_idx = idx;
}

void rastrigin::restoreMove(int idx)
{
    if ((idx != prev_idx) || (!can_rollback))
        throw runtime_error("Rastrigin: Cannot roll back!");
    vars[idx] = prev_x;
    if (outOfBounds)
        outOfBounds = false;
    can_rollback = false;
}

rastrigin::~rastrigin()
{
	if (vars != NULL) {
		delete[] vars;
	}
	vars = NULL;
}
double rastrigin::get_score()
{
    if (outOfBounds)
        return numeric_limits<double>::max();
	double tot = 0;
	int i;
// #pragma omp parallel for private(i) reduction(+:tot)
	for (i = 0; i < dim; i++) {
		tot += vars[i] * vars[i] - 10.0 * cos(2 * M_PI * vars[i]);
	}
	return (10 * dim + tot);
}

double rastrigin::get_param(int idx) const
{
    return vars[idx];
}

void rastrigin::set_param(int idx, double val)
{
	vars[idx] = val;
}

void rastrigin::serialize(void* buf) const
{
    memcpy(buf, vars, sizeof(double) * dim);
}

void rastrigin::deserialize(void const *buf)
{
    memcpy(vars, buf, sizeof(double) * dim);
}

double rastrigin::scramble()
{
    int i;
    for (i = 0; i < dim; ++i) {
        vars[i] = VAR_MAX * (2.0 * rnd.random() - 1.0);
    }
    return get_score();
}
