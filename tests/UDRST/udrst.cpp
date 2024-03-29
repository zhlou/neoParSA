#include "udrst.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <stdexcept>
#include <limits>

#include "unirandom.h"
#include <string.h>
#include <boost/format.hpp>
// #include <omp.h>

using namespace std;
const double udrst::VAR_MAX = 5.12;
const double udrst::VAR_MIN = -5.12;

udrst::udrst(int dimension, unirandom &in_rnd) :
		dim(dimension), rnd(in_rnd)
{
	int i;
	vars = new double[dim];
	for (i = 0; i < dim; i++) {
		vars[i] = VAR_MAX * (rnd.random() *2.0 -1.0);
	}
	//docroot = NULL;
	//section = NULL;
	prev_x = 0; // to make the compiler happy
	prev_idx = -1;
	idx = 0;
	can_rollback = false;
	outOfBounds = false;

}

udrst::udrst(ptree &root, unirandom &in_rnd):
        rnd(in_rnd)
{
    ptree &sec_attr = root.get_child("rastrigin.<xmlattr>");
    dim = sec_attr.get<int>("dim");
    vars = new double[dim];
    for (int i = 0; i < dim; i++) {
        std::string namebuf = (boost::format("x%d") % (i+1)).str();
        boost::optional<double> xval = sec_attr.get_optional<double>(namebuf);
        if (xval)
            vars[i] = *xval;
        else {
            vars[i] = VAR_MAX * (2.0 * rnd.random() - 1.0);
            sec_attr.put(namebuf, boost::format("%g") % vars[i]);
            // for more accuracy, use the following line
            // sec_attr.put(namebuf, vars[i]);
        }
    }
    prev_x = 0; // to make the compiler happy
    prev_idx = -1;
    idx = 0;
    can_rollback = false;
    outOfBounds = false;
}


void udrst::write_section(ptree &root, std::string secname)
{
    ptree section;
    section.put("<xmlattr>.dim", (boost::format("%1%") % dim).str());
	for (int i = 0; i < dim; i++) {
        section.put((boost::format("<xmlattr>.x%1%") % (i+1)).str(),
            (boost::format("%g") % vars[i]).str());
	}
    root.put_child(secname, section);
}

void udrst::print_solution(ostream& o) const
{
	o << "{" << endl;
	for (int i = 0; i < dim; i++) {
		o << vars[i] << endl;
	}
	o << "}" << endl;
}

void udrst::generateMove(int, double theta)
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
    idx = (idx + 1) % dim;
}

void udrst::restoreMove(int)
{
    if (!can_rollback)
        throw runtime_error("Rastrigin: Cannot roll back!");
    vars[prev_idx] = prev_x;
    if (outOfBounds)
        outOfBounds = false;
    can_rollback = false;
}

udrst::~udrst()
{
	if (vars != NULL) {
		delete[] vars;
	}
	vars = NULL;
}
double udrst::get_score()
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

double udrst::get_param(int idx) const
{
    return vars[idx];
}

void udrst::set_param(int idx, double val)
{
	vars[idx] = val;
}

void udrst::serialize(void* buf) const
{
    memcpy(buf, vars, sizeof(double) * dim);
}

void udrst::deserialize(void const *buf)
{
    memcpy(vars, buf, sizeof(double) * dim);
}

double udrst::scramble()
{
    int i;
    for (i = 0; i < dim; ++i) {
        vars[i] = VAR_MAX * (2.0 * rnd.random() - 1.0);
    }
    return get_score();
}
