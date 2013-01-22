#include "tsp.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

void city::print(ostream &o) const
{
    o << "(" << x << "," << y << ")";
}

ostream & operator << (ostream &o, const city &c)
{
    c.print(o);
    return o;
}

double dist(const city &c1, const city &c2)
{
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return sqrt(dx*dx + dy*dy);
}

void tsp::add_city(city the_city)
{
    if (ncities > 0)
        cities.back().next = ncities;
    cities.push_back(the_city);
    can_rollback = false;
    cost_valid = false;
    ncities = cities.size();

}

void tsp::print_array(ostream &o) const
{
    for (vector<city>::const_iterator it = cities.begin(); it != cities.end();
            it++) {
        if (it != cities.begin())
            o << "->";
        o << *it << it->next;
    }
    o << endl;
}

void tsp::print_route(ostream &o) const
{
    if (cities.size() <= 0)
        return;
    
    o << cities[0];
    for (int i = cities[0].next; i != 0; i = cities[i].next) {
        o << "->";
        o << cities[i];
    }
    o << endl;
}

int tsp::get_ncities() const
{
    return cities.size();
}

double tsp::calc_route()
{
    if (cities.size() <=0)
        return 0;
    double cost = 0.0;
    for (int cur=0, i = cities[0].next; i!= 0; cur = i, i = cities[i].next) {
        cost += dist(cities[cur],cities[i]);
    }
    route_cost = cost;
    cost_valid = true;
    return cost;
}

double tsp::get_score()
{
    if (cost_valid)
        return route_cost;
    else
        return calc_route();
}

double tsp::swap(int c1, int c2)
{
    int n1 = cities[c1].next;
    int n2 = cities[c2].next;

    double previous = dist(cities[c1],cities[n1]) 
        + dist(cities[c2], cities[n2]);
    double current = dist(cities[c1],cities[c2]) 
        + dist(cities[n1], cities[n2]);
    double diff = current - previous;

    for (int prev = n1, i = cities[n1].next, next = cities[i].next; i != n2; 
            i = next, next = cities[next].next) {
        assert(i != c1);
        cities[i].next = prev;
        prev = i;
    }
    cities[c1].next = c2;
    cities[n1].next = n2;

    route_cost += diff;
    return diff;

}
double tsp::step(int c1, int c2)
{
    int n = cities.size();
    if (c1 < 0 || c1 >= n || c2 < 0 || c2 >= n)
        return -1;
    if (c1 == c2) {
    	can_rollback = false;
    	return 0;
    }


    r1 = c1;
    r2 = cities[c1].next;
    can_rollback = true;

    return swap(c1, c2);
}

double tsp::roll_back()
{
    if (!can_rollback)
        return 0;

    can_rollback = false;
    return swap(r1, r2);
}

int tsp::getDimension()
{
    return 1;
}

void tsp::generateMove(int, double)
{
    int c1, c2;
    int ncities = cities.size();
    c1 = rand_r(&seed) % ncities;
    do {
        c2 = rand_r(&seed) % ncities;
    } while (c2 == c1);

    step(c1, c2);
}

void tsp::restoreMove(int)
{
    roll_back();
}

tsp::tsp()
{
    can_rollback = false;
    cost_valid = true;
    route_cost = 0.0;
    r1 = r2 = 0; // just to avoid having uninialized variables
    seed = time(NULL);
    ncities = cities.size();
}
