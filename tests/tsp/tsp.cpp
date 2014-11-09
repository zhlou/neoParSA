#include <cassert>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include <libxml/tree.h>

#include "tsp.h"
#include "xmlUtils.h"

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

string tsp::print_route() const
{
	ostringstream convert;
	for (int i = 0; i < ncities; ++i) {
		convert << tour[i] << " -> ";
	}
	convert << tour[0] << endl;
	return convert.str();
}

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
		remove_cost += get_edge(id2,id2p);
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
	can_rollback = true;
    return route_cost;

}

bool tsp::pairLess::operator() (const neighbor_pair& e1, const neighbor_pair& e2) const
{
    return (theTsp.get_edge(e1.from(),e1.to()) < theTsp.get_edge(e2.from(),e2.to()));

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

void tsp::generateMove(int, double theta)
{
    size_t c1, c2;
    c1 = rand_r(&seed) % ncities;
    if (theta < 0)
        theta = - theta;
    c2 = (size_t)theta;
    if ( (c2 >= ncities - 1) || (theta > (double) (ncities -1)))
    	c2 = rand_r(&seed) % (ncities -1);
    //cout << "move (id, theta, neighbor): (" << c1 << ", " << theta << ", " << c2 << ")" << endl;
    c2 = neighbors[c1][c2].to();
    swap(c1,c2);
    // step(c1, c2);
}

void tsp::restoreMove(int)
{
    roll_back();
    can_rollback = false;
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
        for (j = 0; j < i; ++j){
            edge_wt[i][j] = dist(city_list[i],city_list[j]);
        }

    }
    for (i = 0; i < ncities; ++i) {
        for (j = 0; j < ncities; ++j) {
            if (j != i) {
                neighbors[i].push_back(neighbor_pair(i,j));
            }
            sort(neighbors[i].begin(),neighbors[i].end(), pairLess(*this));
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

tsp::tsp(xmlNodePtr docroot)
{
	xmlNodePtr graph = getSectionByName(docroot, "graph");
	if (!graph)
		throw runtime_error(string("Error: fail to find section graph"));

	size_t i = 0, j = 0;
	xmlNodePtr v_iter = graph->children;
	xmlNodePtr e_iter = NULL;
	double wt;
	char * text_buf = NULL;
	for (v_iter = graph->children; v_iter != NULL; v_iter = v_iter->next){
		if (xmlStrcmp(v_iter->name, BAD_CAST "vertex"))
		    continue;
		tour.push_back(i);
		position.push_back(i);
		edge_wt.push_back(vector<double>());
		for (e_iter = v_iter->children; e_iter != NULL; e_iter = e_iter->next) {
		    if (xmlStrcmp(e_iter->name, BAD_CAST "edge"))
		        continue;
		    text_buf = (char *)xmlNodeGetContent(e_iter);
		    sscanf(text_buf,"%lu",&j);
		    xmlFree(text_buf);
		    if (j < i) {
		        text_buf = (char *)xmlGetProp(e_iter,BAD_CAST "cost");
		        sscanf(text_buf,"%lf",&wt);
		        xmlFree(text_buf);
		        edge_wt[i].push_back(wt);
		    }
		}
		++i;
	}
	ncities = tour.size();
	neighbors.resize(ncities);
    for (i = 0; i < ncities; ++i) {
        for (j = 0; j < ncities; ++j) {
            if (j != i) {
                neighbors[i].push_back(neighbor_pair(i,j));
            }
            sort(neighbors[i].begin(),neighbors[i].end(), pairLess(*this));
        }
    }
    can_rollback = false;
    route_cost = calc_tour();
    prev_cost = route_cost;
    r1 = r2 = 0;
    seed = time(NULL);

}

void tsp::save_tsplib_xml(const char *name) const
{
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr root = xmlNewNode(NULL, BAD_CAST
                                 "travellingSalesmanProblemInstance");
    xmlNodePtr graph = xmlNewChild(root, NULL, BAD_CAST "graph",NULL);
    xmlNodePtr vertex_iter = NULL;
    xmlNodePtr edge_iter = NULL;
    char *text_buf = NULL;
    for (size_t i = 0; i < ncities; ++i) {
        vertex_iter = xmlNewChild(graph, NULL, BAD_CAST "vertex", NULL);
        for (size_t j = 0; j < ncities; ++j) {
            if (j != i) {
                asprintf(& text_buf, "%lu",j);
                edge_iter = xmlNewChild(vertex_iter, NULL, BAD_CAST "edge",
                                        BAD_CAST text_buf);
                free(text_buf);
                asprintf(&text_buf, "%.15e",get_edge(i,j));
                xmlNewProp(edge_iter, BAD_CAST"cost", BAD_CAST text_buf);
                free(text_buf);
            }
        }
    }
    xmlDocSetRootElement(doc, root);
    xmlSaveFormatFile(name, doc, 1);
    xmlFreeDoc(doc);
}

void tsp::write_tour(xmlNodePtr xmlroot, const char *tourname)
{
    char text_buf[100];
    xmlNodePtr tourNode = getSectionByName(xmlroot, tourname);
    if (tourNode != NULL) {
        xmlUnlinkNode(tourNode);
        xmlFreeNode(tourNode);
    }
    tourNode = xmlNewChild(xmlroot, NULL, BAD_CAST tourname, NULL);
    for (size_t i = 0; i < ncities; ++i) {
        sprintf(text_buf,"%lu",tour[i]);
        xmlNewChild(tourNode, NULL, BAD_CAST "vertex", BAD_CAST text_buf);
    }
}

double tsp::read_tour(const xmlNodePtr xmlroot, const char *tourname)
{
    xmlNodePtr xmltour = getSectionByName(xmlroot, tourname);
    if (NULL == xmltour) {
        throw runtime_error(string("Error: fail to find tour section"));
    }
    vector<size_t> tmp_tour(ncities,0);
    vector<size_t> tmp_position(ncities,-1);
    xmlChar *content = NULL;
    size_t i = 0, j;
    for (xmlNodePtr vNode = xmltour->children; vNode != NULL; vNode = vNode->next) {
        if (xmlStrcmp(vNode->name, BAD_CAST"vertex"))
            continue;
        if (i >= ncities) // ignore
            break;
        content = xmlNodeGetContent(vNode);
        if (content == NULL)
            throw runtime_error(string("Error: tour vertex with no ID"));
        tmp_tour[i] = atoi((char *)content);
        xmlFree(content);
        ++i;
    }
    if (i < ncities)
        throw runtime_error(string("Error: not enough vertices in tour"));
    for (i = 0; i < ncities; ++i) {
        j = tmp_tour[i];
        if (tmp_position[j] != -1)
            throw runtime_error(string("invalid tour"));
        tmp_position[j] = i;
    }
    copy(tmp_tour.begin(),tmp_tour.end(),tour.begin());
    copy(tmp_position.begin(),tmp_position.end(),position.begin());

    can_rollback = false;
    return (route_cost = calc_tour());
}

void tsp::serialize(void *buf) const
{
    double *score = (double *)buf;
    *score = route_cost;
	size_t *buf_tour = (size_t*)((double*)buf+1);
	for (size_t i = 0; i < ncities; ++i) {
		buf_tour[i] = tour[i];
	}

}

void tsp::deserialize(const void *buf)
{
    route_cost = *(double *)buf;
    size_t *buf_tour = (size_t*)((double*)buf+1);
	size_t i;
	for (i = 0; i < ncities; ++i) {
		tour[i] = buf_tour[i];
		position[tour[i]] = i;
	}
	can_rollback = false;
}
