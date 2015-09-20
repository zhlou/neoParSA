/* 
 * File:   vectorUtils.h
 * Author: zhlou
 *
 * Created on September 19, 2015, 8:23 PM
 */

#ifndef VECTORUTILS_H
#define	VECTORUTILS_H

#include <vector>
#include <fstream>

template<typename T>
void readVector(std::vector<T> &dat, std::ifstream &is)
{
    typename std::vector<T>::size_type sz;
    is.read((char *)&sz, sizeof(sz));
    dat.resize(sz);
    is.read((char*)&dat[0], sz*sizeof(T));
}


template<typename T>
void readVector(std::vector<T> &dat, const char* filename)
{
    std::ifstream is(filename, std::ios::binary);
    readVector(dat, is);
    is.close();
}

template<typename T>
void writeVector(const std::vector<T> &dat, std::ofstream &os)
{
    typename std::vector<T>::size_type sz = dat.size();
    os.write((const char*)&sz, sizeof(sz));
    os.write((const char*)&dat[0], sz*sizeof(T));
}

template<typename T>
void writeVector(const std::vector<T> &dat, const char* filename)
{
    std::ofstream os(filename, std::ios::binary);
    writeVector(dat, os);
    os.close();
}

void readDoubleVectorFromText(std::vector<double> &dat, const char* filename)
{
    std::ifstream is(filename);
    double val;
    while (is >> val)
        dat.push_back(val);
    is.close();
}

#endif	/* VECTORUTILS_H */

