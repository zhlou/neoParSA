
/*
 * xmlUtils.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: zhlou
 */

#include "xmlUtils.h"
#include <string>
#include <stdexcept>
#include <cstdlib>
using namespace std;

// return double value from node "section" with property name "name"
// throws runtime_error when not found
double getPropDouble(xmlNode* section, const char* name)
{
    if (section == NULL)
        throw runtime_error(string("Error: NULL section"));
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (xmlChar *)name);
    if (prop == NULL)
        throw runtime_error(string("Error: fail to find property ")
                +name+" in xml section "+(char *)(section->name));
    double val = strtod((char *)prop, NULL);
    xmlFree(prop);
    return val;
}

// return long value from node "section" with property name "name"
// throws runtime_error when not found
long getPropLong(xmlNode* section, const char* name)
{
    if (section == NULL)
        throw runtime_error(string("Error: NULL section"));
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (xmlChar *)name);
    if (prop == NULL)
        throw runtime_error(string("Error: fail to find property ")
                +name+" in xml section "+(char *)(section->name));
    long val = strtol((char *)prop, NULL, 0);
    xmlFree(prop);
    return val;
}

// return integer value from node "section" with property name "name"
// throws runtime_error when not found
int getPropInt(xmlNode* section, const char* name)
{
    if (section == NULL)
            throw runtime_error(string("Error: NULL section"));
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (const xmlChar *)name);
    if (prop == NULL)
        throw runtime_error(string("Error: fail to find property ")
                +name+" in xml section "+(char *)(section->name));
    int val = atoi((char *)prop);
    xmlFree(prop);
    return val;
}

// return node with name "name" or NULL if not found
xmlNode* getSectionByName(xmlNode* root, const char* name)
{
    xmlNode *section = root->children;
    while (section != NULL) {
        if(!xmlStrcmp(section->name,(const xmlChar *)name))
            break;
        section = section->next;
    }
    return section;
}
