/*
 * flyUtils.h
 *
 *  Created on: Feb 4, 2013
 *      Author: zhlou
 */

#ifndef FLYUTILS_H_
#define FLYUTILS_H_

const int MAX_RECORD = 256;
const int MAX_ARGS = 25;

void error(const char *format, ...);
void warning(const char *format, ...);
FILE *FindSection(FILE *fp, const char *section);


#endif /* FLYUTILS_H_ */
