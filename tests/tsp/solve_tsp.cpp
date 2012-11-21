#include <iostream>
#include <fstream>
#include "annealer.h"
#include "movable.h"
#include "tsp.h"
#include "tsp_problem.h"

using namespace std;

int main(int argc, char **argv)
{
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
    annealer tsp_anneal(&test_problem);
    cout << "The finial energy is "<<tsp_anneal.loop()<< endl;
    cout << "The energy cached is "<< test_tsp.cost() << endl;
    cout << "The real energy is " << test_tsp.calc_route()<<endl;

    test_tsp.print_route(cout);
    return 0;

}
