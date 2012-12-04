/*
 * test_rst.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: zhlou
 */

#include "rastrigin.h"
#include <libxml/parser.h>
#include <iostream>

int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr << "Missing input files" << endl;
	}
	char *docname = argv[1];
	xmlDoc *doc = xmlParseFile(docname);
	xmlNode *root = xmlDocGetRootElement(doc);
	rastrigin rst(root);
	rst.write_section((xmlChar *)"output");
	xmlSaveFile(docname, doc);
	xmlFree(doc);
	return 0;
}


