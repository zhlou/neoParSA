/*
 * xmlUtils.h
 *
 *  Created on: Dec 14, 2012
 *      Author: zhlou
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <libxml/tree.h>
double getPropDouble(xmlNode *section, const char *name);
int getPropInt(xmlNode *section, const char *name);
long getPropLong(xmlNode *section, const char *name);
xmlNode *getSectionByName(xmlNode *root, const char *name);


#endif /* UTILS_H_ */
