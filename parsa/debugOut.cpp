/*
 * debugOut.cpp
 *
 *  Created on: Jan 29, 2013
 *      Author: zhlou
 */


#include "debugOut.h"


//
//ignore& ignore::operator <<(StandardEndLine char_traits)
//{
//    return *this;
//}

ignore debugIGNORE::debugOut;
ostream &debugSTD::debugOut = cout;
ofstream debugNULL::debugOut("/dev/null");
