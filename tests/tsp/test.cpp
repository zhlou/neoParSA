#include "tsp.h"
#include <iostream>

using namespace std;

int main (int argc, char **argv)
{
    vector<city> cities;
    cities.push_back(city(0,1));
    cities.push_back(city(0,2));
    cities.push_back(city(0,3));
    cities.push_back(city(0,4));
    vector<city>::iterator it;
    for (it = cities.begin(); it != cities.end(); it++) {
        cout << *it << endl;
    }

    tsp the_problem;
    the_problem.add_city(city(1,1));
    the_problem.add_city(city(1,2));
    the_problem.add_city(city(1,3));
    the_problem.add_city(city(1,4));
    the_problem.add_city(city(1,5));
    the_problem.print_array(cout);
    the_problem.print_route(cout);
    cout << "The cost of the route is" << the_problem.calc_route() << endl;
    cout << the_problem.step(1,3) << endl;
    the_problem.print_route(cout);
    cout << "The new cost is" << the_problem.calc_route() << endl;
    return 0;
}


