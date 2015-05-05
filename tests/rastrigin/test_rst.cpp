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
#include <libxml/parser.h>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr << "Missing input files" << endl;
		return 1;
	}
	unirandom rnd;
	char *docname = argv[1];
	xmlDoc *doc = xmlParseFile(docname);
	xmlNode *root = xmlDocGetRootElement(doc);
	rastrigin rst(root, rnd);
	cout << "The score is " << rst.get_score() << endl;
	double nScore = rst.scramble();
	cout << "New score is " << nScore << endl;
	rst.print_solution(cout);
	rst.write_section((xmlChar *)"rastrigin");
	expHold::Param scheduleParam(root);
	tempCount::Param tmpCntParam(root);
	annealer<rastrigin, expHold, tempCount, feedbackMove>
	    rst_anneal(rst, rnd, scheduleParam, tmpCntParam, root);
	std::string steplogName(argv[1]);
	steplogName.replace(steplogName.end()-4,steplogName.end(),".steplog");
	rst_anneal.setStepLog(file, steplogName.c_str());
	rst_anneal.loop();
	xmlSaveFormatFile(docname, doc,1);
	xmlFree(doc);
	return 0;
}


