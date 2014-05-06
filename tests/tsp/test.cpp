#include "tsp.h"
#include <fstream>
#include <iostream>

using namespace std;

int main (int argc, char **argv)
{
    tsp test_tsp;
//    ifstream ifile("input.dat");
//    if (!ifile)
//        return -1;
//    double x,y;
//    while (!ifile.eof())
//    {
//        ifile >> x >> y;
//        test_tsp.add_city(city(x,y));
//    }
    test_tsp.add_city(city(0,0));//0
    test_tsp.add_city(city(1,0));//1
    test_tsp.add_city(city(2,0));//2
    test_tsp.add_city(city(3,0));//3
    test_tsp.add_city(city(3,1));//4
    test_tsp.add_city(city(2,1));//5
    test_tsp.add_city(city(1,1));//6
    test_tsp.add_city(city(0,1));//7

    cout.precision(16);
    test_tsp.print_route(cout);
    cout << "The initial cost is "<<test_tsp.get_score() << endl;
    cout << "The real initial cost is "<<test_tsp.calc_route() << endl;
    double cost = test_tsp.get_score();
    cost += test_tsp.step(0, 7);
    cout << "The cost after step is "<< cost << endl;
    cout << "and it should be "<< test_tsp.get_score() << endl;
    test_tsp.print_route(cout);
    test_tsp.roll_back();
    cout << "The cost after rollback is "<<test_tsp.get_score() << endl;
    cout << "The real cost is "<<test_tsp.calc_route() << endl;
    test_tsp.print_route(cout);

    return 0;
}



