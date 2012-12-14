/*
 * utils.h
 *
 *  Created on: Dec 14, 2012
 *      Author: zhlou
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <libxml/tree.h>
double getPropDouble(xmlNode *section, char *name);
int getPropInt(xmlNode *section, char *name);
xmlNode *getSectionByName(xmlNode *root, char *name);


#endif /* UTILS_H_ */
