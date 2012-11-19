#ifndef TSP_H
#define TSP_H

#include<vector>
#include<list>
#include<ostream>

using namespace std;



class city {
    private:
        double x,y;
    public:
        int next;
        city(double _x, double _y) : x(_x), y(_y), next(0){}
        friend double dist(const city&, const city&);
        void print(ostream &) const;
};
ostream & operator << (ostream &, const city &);
double dist(const city &c1, const city &c2);

class tsp {
    private:
        vector<city> cities;
        bool can_rollback;
        bool cost_valid;
        double route_cost;
        int r1, r2;
        double swap(int c1, int c2);
    public:
        void print_array(ostream &) const;
        void print_route(ostream &) const;
        void add_city(city);
        int get_ncities() const;
        double calc_route();
        double cost();
        double step(int c1, int c2);
        double roll_back();
        tsp();
};

#endif
