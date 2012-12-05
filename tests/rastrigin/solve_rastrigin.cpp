#include <iostream>
#include <libxml/parser.h>
#include "rastrigin.h"
#include "rastrigin_problem.h"
#include "unirandom.h"
#include "annealer.h"
using namespace std;
int main(int argc, char **argv)
{
	if (argc <= 1) {
		cerr<< "Missing input files" << endl;
		return 1;
	}
	char *docname = argv[1];
	xmlDoc *doc = xmlParseFile(docname);
	xmlNode *docroot = xmlDocGetRootElement(doc);
	if (doc == NULL) {
		cerr << "Input incorrect" << endl;
		return 2;
	}
	unirandom rnd;
	rastrigin rst(3, rnd);
	rastrigin_problem rst_problem(&rst);
	annealer rst_anneal(&rst_problem, docroot);
	cout << "The initial state is " << endl;
	rst.print_solution(cout);
	cout << "The fininal energy is " << rst_anneal.loop() << endl;
	cout << "The solution is " << endl;
	rst.print_solution(cout);

	return 0;
}
