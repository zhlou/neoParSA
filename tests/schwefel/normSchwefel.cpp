/*
 * normSchwefel.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: zhlou
 */

#include "normSchwefel.h"

#include <limits>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <libxml/tree.h>

#include "unirandom.h"
#include "xmlUtils.h"

const double normSchwefel::VAR_MAX = 512.0;
const double normSchwefel::VAR_MIN = -512.0;

normSchwefel::normSchwefel(int dimension, unirandom& in_rnd) :
        dim(dimension), rnd(in_rnd)
{
    int i;
    vars = new double[dim];
    for (i = 0; i < dim; ++i) {
        vars[i] = VAR_MAX * (rnd.random() * 2.0 - 1.0);
    }
    docroot = NULL;
    section = NULL;
    prev_x = 0;
    prev_idx = -1;
    idx = 0;
    can_rollback = false;
    outOfBounds = false;
}

normSchwefel::normSchwefel(xmlNode* root, unirandom& in_rnd,
                           const char* secName) :
        docroot(root), rnd(in_rnd)
{
    section = getSectionByName(docroot, secName);
    if (section == NULL)
        throw 1;
    dim = getPropInt(section, "dim");
    vars = new double[dim];
    char *namebuf = NULL;
    for (int i = 0; i < dim; i++) {
        asprintf(&namebuf, "x%d", i);
        try {
            vars[i] = getPropDouble(section, namebuf);
        } catch (std::exception &e) {
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

    prev_x = 0;
    prev_idx = -1;
    idx = 0;
    can_rollback = false;
    outOfBounds = false;
}

normSchwefel::~normSchwefel()
{
    if (vars != NULL) {
        delete[] vars;
    }
    vars = NULL;
}

double normSchwefel::getParam(int idx) const
{
    return vars[idx];
}

void normSchwefel::setParam(int idx, double val)
{
    vars[idx] = val;
}

double normSchwefel::get_score()
{
    if (outOfBounds)
        return std::numeric_limits<double>::max();
    double tot = 0.0;
    double x;
    for (int i = 0; i < dim; ++i) {
        x = vars[i];
        tot -= x * std::sin(std::sqrt(std::abs(x)));
    }
    return (tot / (double)dim);
}

void normSchwefel::generateMove(int, double theta)
{
    double x = prev_x = vars[idx];
    x += theta;
    if (x > VAR_MAX || x < VAR_MIN) {
        x = prev_x;
        outOfBounds = true;
    }
    vars[idx] = x;
    can_rollback = true;
    prev_idx = idx;
    ++ idx;
    idx = (dim == idx) ? 0 : idx;
}

void normSchwefel::restoreMove(int)
{
    if (!can_rollback)
        throw std::runtime_error("normSchwefel: Cannot roll back!");
    vars[prev_idx] = prev_x;
    outOfBounds = false;
    can_rollback = false;
}

void normSchwefel::printSolution(std::ostream& o) const
{
    o << "{\n";
    for (int i = 0; i < dim; ++i) {
        o << vars[i] << '\n';
    }
    o << "}" << std::endl;
}

void normSchwefel::writeSection(xmlChar* secname)
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
        sprintf(namebuf, "x%d", i);
        sprintf(valbuf, "%g", vars[i]);
        xmlNewProp(node, (xmlChar *) namebuf, (xmlChar *) valbuf);
    }
    delete[] valbuf;
    delete[] namebuf;
}

void normSchwefel::serialize(void* buf) const
{
    std::memcpy(buf, vars, sizeof(double) * dim);
}

void normSchwefel::deserialize(const void* buf)
{
    std::memcpy(vars, buf, sizeof(double) * dim);
}

double normSchwefel::scramble()
{
    for (int i = 0; i < dim; ++i) {
        vars[i] = VAR_MAX * (2.0 * rnd.random() - 1.0);
    }
    return get_score();
}
