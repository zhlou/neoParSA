#include <iostream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
using boost::property_tree::ptree;
#include "rastrigin.h"
#include "annealer.h"
//#include "simpleAnnealer.h"
#include "move/feedbackMove.h"
//#include "rastrigin_problem.h"
#include "unirandom.h"
#include "lam.h"
#include "criCount.h"
using namespace std;
int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr<< "Missing input files" << endl;
		return 1;
	}
	char *docname = argv[1];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &docroot = pt.begin()->second;
	unirandom rnd;
	rastrigin rst(docroot, rnd);
	//feedbackMove<rastrigin, debugSTD> rst_problem(rst, rnd, docroot);
	//lam schedule(docroot);
	lam::Param scheduleParam(docroot);
	criCount::Param criCntParam(docroot);
	annealer<rastrigin, lam, criCount, feedbackMove>
	        rst_anneal(rst, rnd, scheduleParam, criCntParam, docroot);
	//rastrigin_problem rst_problem(&rst, docroot);
	//lam rst_anneal(&rst_problem, docroot);
	cout << "The initial state is " << endl;
	rst.print_solution(cout);
	cout << "The fininal energy is " << rst_anneal.loop() << endl;
	cout << "The solution is " << endl;
	rst.print_solution(cout);
	rst_anneal.writeResult(docroot);
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(docname, pt, std::locale(), settings);

	return 0;
}
