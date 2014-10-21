#include <fstream>
#include <iostream>
#include <vector>

#include <libxml/parser.h>
#include "tsp.h"

using namespace std;

int main (int argc, char **argv)
{
    char xmlname[] = "test8.xml";
    xmlDocPtr doc = xmlParseFile(xmlname);
    xmlNodePtr root = xmlDocGetRootElement(doc);
    tsp test_tsp(root);
    /*
    vector<city> clist;
    clist.push_back(city(0,1));
    clist.push_back(city(0,0));
    clist.push_back(city(1,1));
    clist.push_back(city(1,0));
    clist.push_back(city(2,1));
    clist.push_back(city(2,0));
    clist.push_back(city(3,1));
    clist.push_back(city(3,0));

    tsp test_tsp(clist);
    */
    test_tsp.save_tsplib_xml("test8.xml");

    cout.precision(16);

    cout << "The initial cost is " << test_tsp.get_score() << endl;
    cout << "The initail tour is:" << endl << test_tsp.print_route();
    test_tsp.swap(2,3);
    cout << endl;
    cout << "After swap, cost is " << test_tsp.get_score() << endl;
    cout << "            real cost is " << test_tsp.calc_tour() << endl;
    cout << "            tour is " << endl << test_tsp.print_route();
//    ifstream ifile("input.dat");
//    if (!ifile)
//        return -1;
//    double x,y;
//    while (!ifile.eof())
//    {
//        ifile >> x >> y;
//        test_tsp.add_city(city(x,y));
//    }
    /*
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
    double cost;
    for (int i = 0; i < 8; ++i){
        for (int j = 0; j < 8; ++j){
            if (j != i) {
                cout << "swapping " << i << " and " << j << endl;
                cout << "The initial cost is "<<test_tsp.get_score() << endl;
                cout << "The real initial cost is "<<test_tsp.calc_route() << endl;
                cost = test_tsp.get_score();
                cost += test_tsp.step(i, j);
                cout << "The cost after step is "<< cost << endl;
                cout << "and it should be "<< test_tsp.get_score() << endl;
                test_tsp.print_route(cout);
                test_tsp.roll_back();
                cout << "The cost after rollback is "<<test_tsp.get_score() << endl;
                cout << "The real cost is "<<test_tsp.calc_route() << endl;
                test_tsp.print_route(cout);
            }
        }
    }

*/
    return 0;
}



