#ifndef TSP_H
#define TSP_H

#include <vector>
#include <list>
#include <ostream>
#include <string>

#include <libxml/tree.h>

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
    // This is the functor used in sorting the neighbors.
    struct pairLess
    {
        tsp& theTsp;
        pairLess(tsp& theTsp):theTsp(theTsp){}
        bool operator () (const neighbor_pair& e1, const neighbor_pair& e2) const;
    };
    /*
     * Make sure that tour[position[i]] == i
     */
    vector<size_t> tour;
    vector<size_t> position;
    vector<vector <double> > edge_wt;
    vector<vector <neighbor_pair > > neighbors;

    bool can_rollback;
    double prev_cost;
    double route_cost;
    int r1, r2;
    unsigned seed;
    size_t ncities;

    double get_edge(size_t i, size_t j) const {return i > j ? edge_wt[i][j] : edge_wt[j][i];}
    size_t prev(size_t i){ return (i>0)?(i-1):(ncities-1);}
    size_t next(size_t i){ return (i<ncities-1)?(i+1):0;}

    tsp(const tsp &tsp); // = delete; disable copy constructor until it gets implemented
public:
    tsp();
    tsp(vector<city>& city_list);
    tsp(xmlNodePtr docroot);
    string print_route() const;
    size_t get_ncities() const {return ncities;}
    double swap(size_t c1, size_t c2);
    double calc_tour() const;
    int getDimension();
    void generateMove(int, double);
    void restoreMove(int);
    double get_score() const {return route_cost;}
    double roll_back();
    void save_tsplib_xml(const char* name) const;
    void write_tour(xmlNodePtr xmlroot, const char *tourname);
    double read_tour(const xmlNodePtr xmlroot, const char *tourname);
    size_t getStateSize() const {return sizeof(size_t)*ncities+sizeof(double);}
    void serialize(void *buf) const;
    void deserialize(void const *buf);
};

#endif
