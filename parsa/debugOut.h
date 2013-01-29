/*
 * debugOut.h
 *
 *  Created on: Jan 29, 2013
 *      Author: zhlou
 */

#ifndef DEBUGOUT_H_
#define DEBUGOUT_H_
#include <iostream>
#include <fstream>
using namespace std;
class ignore
{
public:
    template<class T>
    ignore& operator <<(T const &);


    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    ignore& operator<<(StandardEndLine) {return *this;}

};
//template<class T>
//ignore& operator << (ignore &ig, T const &);
class debugIGNORE
{
protected:
    static ignore debugOut;
};

class debugSTD
{
protected:
    static ostream &debugOut;
};

class debugNULL
{
protected:
    static ofstream debugOut;
};

template<class T>
inline ignore& ignore::operator <<(T const &)
{
    return *this;
}
#endif /* DEBUGOUT_H_ */
