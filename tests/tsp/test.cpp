#include "tsp.h"
#include <fstream>
#include <iostream>

using namespace std;

int main (int argc, char **argv)
{
    tsp test_tsp;
    ifstream ifile("input.dat");
    if (!ifile)
        return -1;
    double x,y;
    while (!ifile.eof())
    {
        ifile >> x >> y;
        test_tsp.add_city(city(x,y));
    }
    cout.precision(16);
    test_tsp.print_route(cout);
    cout << "The initial cost is "<<test_tsp.cost() << endl;
    cout << "The real initial cost is "<<test_tsp.calc_route() << endl;
    double cost = test_tsp.cost();
    cost += test_tsp.step(3, 7);
    cout << "The cost after step is "<< cost << endl;
    cout << "and it should be "<< test_tsp.cost() << endl;
    test_tsp.roll_back();
    cout << "The cost after rollback is "<<test_tsp.cost() << endl;
    cout << "The real cost is "<<test_tsp.calc_route() << endl;


    return 0;
}


