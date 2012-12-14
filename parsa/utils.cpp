
/*
 * utils.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: zhlou
 */

#include "utils.h"
#include <stdexcept>
#include <cstdlib>
using namespace std;


double getPropDouble(xmlNode* section, char* name)
{
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (xmlChar *)name);
    if (prop == NULL)
        throw runtime_error(string("Error: fail to find property ")
                +name+" in xml section "+section->name);
    double val = strtod((char *)prop, NULL);
    xmlFree(prop);
    return val;
}

int getPropInt(xmlNode* section, char* name)
{
    xmlChar *prop = NULL;
    prop = xmlGetProp(section, (xmlChar *)name);
    if (prop == NULL)
        throw runtime_error(string("Error: fail to find property ")
                +name+" in xml section "+section->name);
    int val = atoi((char *)prop);
    xmlFree(prop);
    return val;
}

xmlNode* getSectionByName(xmlNode* root, char* name)
{
    xmlNode *section = root->children;
    while (section != NULL) {
        if(!xmlStrcmp(section->name,(xmlChar *)name))
            break;
        section = section->next;
    }
    return section;
}
