#include <iostream>
#include <fstream>
#include <libxml/parser.h>
#include "simpleAnnealer.h"
#include "movable.h"
#include "tsp.h"
#include "tsp_problem.h"

using namespace std;

int main(int argc, char **argv)
{
	const char *xmlfile = "tsp.xml";
	xmlDoc *doc = xmlParseFile(xmlfile);
	xmlNode *root = xmlDocGetRootElement(doc);
    tsp test_tsp;
    ifstream ifile("input.dat");
    if (!ifile)
        return -1;
    double x,y;
    while(!ifile.eof())
    {
        ifile >> x >> y;
        test_tsp.add_city(city(x,y));
    }
    tsp_problem test_problem(&test_tsp);
    test_tsp.print_route(cout);
    simpleAnnealer tsp_anneal(&test_problem, root);
    cout << "The finial energy is "<<tsp_anneal.loop()<< endl;
    cout << "The energy cached is "<< test_tsp.cost() << endl;
    cout << "The real energy is " << test_tsp.calc_route()<<endl;

    test_tsp.print_route(cout);
    return 0;

}
