/*
 * test_rst.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: zhlou
 */

#include "rastrigin.h"
#include "unirandom.h"
#include <libxml/parser.h>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr << "Missing input files" << endl;
	}
	unirandom rnd;
	char *docname = argv[1];
	xmlDoc *doc = xmlParseFile(docname);
	xmlNode *root = xmlDocGetRootElement(doc);
	rastrigin rst(root, rnd);
	cout << "The score is " << rst.get_score() << endl;
	rst.write_section((xmlChar *)"output");
	xmlSaveFormatFile(docname, doc,1);
	xmlFree(doc);
	return 0;
}


