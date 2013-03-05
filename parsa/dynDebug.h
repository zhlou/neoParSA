// This is the base class for dynamic debug output. The purpose of having this
// class is to unify the debug output to stdout, a file or no output at all, in
// iostream like operation.
#ifndef DYNDEBUG_H_
#define DYNDEBUG_H_

#include <ostream>
#include <fstream>

class dynDebug
{
private:
    bool isDebug;
    std::ostream *out;
public:
    dynDebug(int i = 0);
    dynDebug(const char *outname) : isDebug(true), out(new std::ofstream(outname)) {}
    ~dynDebug(){if (out && out != &std::cout && out != &std::cerr) delete out;}
    template<class T>
    dynDebug& operator <<(T const &);



    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    dynDebug& operator<<(StandardEndLine);
};

template<class T>
dynDebug& dynDebug::operator<<(const T& t)
{
    if (isDebug)
        *out << t;
    return *this;
}

#endif
