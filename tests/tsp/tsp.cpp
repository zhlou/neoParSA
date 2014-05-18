#include "tsp.h"
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>

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

/*
void tsp::add_city(city the_city)
{
    if (ncities > 0)
        cities.back().next = ncities;
    cities.push_back(the_city);
    can_rollback = false;
    ncities = cities.size();
    cities.back().next = 0;

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
*/
double tsp::calc_tour() const
{
    double cost = .0;
    if (ncities < 2)
        return cost;
    size_t i,c1,c2;
    for (i = 1; i < ncities; ++i) {
        cost += get_edge(tour[i-1],tour[i]);
    }
    cost += get_edge(tour[ncities-1],tour[0]);

    return cost;
}

double tsp::swap(size_t id1, size_t id2)
{
	size_t p1 = position[id1];
	size_t p2 = position[id2];
	size_t id1p = tour[prev(p1)];
	size_t id1n = tour[next(p1)];
	size_t id2p = tour[prev(p2)];
	size_t id2n = tour[next(p2)];
	double add_cost = 0.0;
	double remove_cost = 0.0;
	prev_cost = route_cost;
	r1 = id1;
	r2 = id2;
	if (p1 == prev(p2)) {
		assert(id1n == id2);
		assert(id1 == id2p);
		remove_cost += get_edge(id1,id1p);
		remove_cost += get_edge(id2,id2n);
		add_cost += get_edge(id1,id2n);
		add_cost += get_edge(id2,id1p);
	} else if (p2 == prev(p1)) {
		assert(id2n == id1);
		assert(id2 == id1p);
		remove_cost += get_edge(id1,id1n);
		remove_cost += get_edge(id2,id1p);
		add_cost += get_edge(id1,id2p);
		add_cost += get_edge(id2,id1n);
	} else {
		assert(id1 != id2p);
		assert(id1 != id2n);
		assert(id2 != id1p);
		assert(id2 != id1n);
		remove_cost += get_edge(id1, id1p);
		remove_cost += get_edge(id1, id1n);
		remove_cost += get_edge(id2, id2p);
		remove_cost += get_edge(id2, id2n);
		add_cost += get_edge(id1, id2p);
		add_cost += get_edge(id1, id2n);
		add_cost += get_edge(id2, id1p);
		add_cost += get_edge(id2, id1n);
	}

	tour[p1] = id2;
	tour[p2] = id1;
	position[id1] = p2;
	position[id2] = p1;
	route_cost += add_cost - remove_cost;
    return route_cost;

}

bool tsp::pairLess::operator() (const neighbor_pair& e1, const neighbor_pair& e2) const
{
    return (theTsp.get_edge(e1.from(),e1.to()) < theTsp.get_edge(e2.from(),e2.to()));

}

/*
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
*/
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

void tsp::generateMove(int, double theta)
{
    size_t c1, c2;
    c1 = rand_r(&seed) % ncities;
    c2 = c1 + (size_t)theta + 1;
    if (c2 >= ncities - 1)
    	c2 = rand_r(&seed) % (ncities -1);
    c2 = neighbors[c1][c2].to();
    swap(c1,c2);
    // step(c1, c2);
}

void tsp::restoreMove(int)
{
    roll_back();
}

tsp::tsp(vector<city>& city_list)
{
    size_t i,j;
    ncities = city_list.size();
    double dx, dy;
    tour.resize(ncities);
    position.resize(ncities);
    edge_wt.resize(ncities);
    neighbors.resize(ncities);
    for (i=0; i < ncities; ++i){
        tour[i] = i;
        position[i] = i;
        edge_wt[i].resize(i);
        //edge_wt[i].resize(i,0);
        for (j = 0; j < i; ++j){
            edge_wt[i][j] = dist(city_list[i],city_list[j]);
        }

    }
    for (i = 0; i < ncities; ++i) {
        for (j = 0; j < ncities; ++j) {
            if (j != i) {
                neighbors[i].push_back(neighbor_pair(i,j));
                sort(neighbors[i].begin(),neighbors[i].end(), pairLess(*this));
            }
        }
    }
    can_rollback = false;
    route_cost = calc_tour();
    prev_cost = route_cost;
    r1 = r2 = 0;
    seed = time(NULL);
}

tsp::tsp()
{
    can_rollback = false;
    route_cost = 0.0;
    prev_cost = route_cost;
    r1 = r2 = 0; // just to avoid having uninialized variables
    seed = time(NULL);
    ncities = 0;
}
