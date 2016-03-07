/*
 * test_rst.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: zhlou
 */

#include "rastrigin.h"
#include "unirandom.h"
#include "annealer.h"
#include "tempCount.h"
#include "expHold.h"
#include "move/feedbackMove.h"
#include <string>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr << "Missing input files" << endl;
		return 1;
	}
	unirandom rnd;
	char *docname = argv[1];
    ptree pt;
    read_xml(docname, pt, boost::property_tree::xml_parser::trim_whitespace);
    ptree &root = pt.begin()->second;
	rastrigin rst(root, rnd);
	cout << "The score is " << rst.get_score() << endl;
	double nScore = rst.scramble();
	cout << "New score is " << nScore << endl;
	rst.print_solution(cout);
	rst.write_section(root, "rastrigin");
	expHold::Param scheduleParam(root);
	tempCount::Param tmpCntParam(root);
	annealer<rastrigin, expHold, tempCount, feedbackMove>
	    rst_anneal(rst, rnd, scheduleParam, tmpCntParam, root);
	std::string steplogName(argv[1]);
	steplogName.replace(steplogName.end()-4,steplogName.end(),".steplog");
	rst_anneal.setStepLog(file, steplogName.c_str());
	rst_anneal.loop();
    boost::property_tree::xml_writer_settings<std::string> settings(' ', 2);
    write_xml(docname, pt, std::locale(), settings);
	return 0;
}
