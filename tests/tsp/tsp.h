#ifndef TSP_H
#define TSP_H

#include<vector>
#include<list>
#include<ostream>

using namespace std;

class city {
private:
    double x, y;
public:
    //int next;
    city(double _x, double _y) :
            x(_x), y(_y)
    {
    }
    friend double dist(const city&, const city&);
    void print(ostream &) const;
};
ostream & operator <<(ostream &, const city &);
double dist(const city &c1, const city &c2);

class tsp {
private:
    struct neighbor_pair
    {
        size_t pair[2];
        neighbor_pair(size_t from, size_t to){pair[0]=from; pair[1]=to;}
        size_t from() const {return pair[0];}
        size_t& from() {return pair[0];}
        size_t to() const {return pair[1];}
        size_t& to() {return pair[1];}
    };
    struct pairLess
    {
        tsp& theTsp;
        pairLess(tsp theTsp):theTsp(theTsp){}
        bool operator () (const neighbor_pair& e1, const neighbor_pair& e2) const;
    };
    vector<size_t> tour;
    vector<vector <double> > edge_wt;
    vector<vector <neighbor_pair > > neighbors;
    //vector<city> cities;
    bool can_rollback;
    double prev_cost;
    //mutable bool cost_valid;
    double route_cost;
    int r1, r2;
    double swap(size_t c1, size_t c2);
    unsigned seed;
    size_t ncities;
    double get_edge(size_t i, size_t j) const {return i > j ? edge_wt[i][j] : edge_wt[j][i];}
    //bool less(const neighbor_pair& e1, const neighbor_pair& e2) const;
public:
    tsp(vector<city>& city_list);
    tsp();
    //void print_array(ostream &) const;
    //void print_route(ostream &) const;
    // void add_city(city);
    size_t get_ncities() const {return ncities;}
    double calc_tour() const;
    int getDimension();
    void generateMove(int, double);
    void restoreMove(int);
    double get_score() const {return route_cost;}
    //double step(int c1, int c2);
    double roll_back();
};

#endif
