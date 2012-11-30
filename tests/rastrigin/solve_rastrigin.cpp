#include <iostream>
#include <libxml/parser.h>
#include "rastrigin.h"
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
	xmlNode *root = xmlDocGetRootElement(doc);
	if (doc == NULL) {
		cerr << "Input incorrect" << endl;
		return 2;
	}
	rastrigin rst_problem(3);
	annealer rst_anneal(&rst_problem, root);
	cout << "The initial state is " << endl;
	rst_problem.print_solution(cout);
	cout << "The fininal energy is " << rst_anneal.loop() << endl;
	cout << "The solution is " << endl;
	rst_problem.print_solution(cout);

	return 0;
}
