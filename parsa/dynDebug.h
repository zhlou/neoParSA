// This is the base class for dynamic debug output. The purpose of having this
// class is to unify the debug output to stdout, a file or no output at all, in
// iostream like operation.
#ifndef DYNDEBUG_H_
#define DYNDEBUG_H_

#include <ostream>
#include <fstream>
enum debugStatus {ignore, out, error, file};
class dynDebug
{
private:
    std::ostream *streamout;
    debugStatus status;
public:
    dynDebug(debugStatus st=ignore, const char *outname=NULL);
    ~dynDebug(){if (status == file) delete streamout;}
    template<class T>
    dynDebug& operator <<(T const &);
    void setDebug(debugStatus st, const char *outname=NULL);



    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    dynDebug& operator<<(StandardEndLine);

protected:
    void __setDebug(const char* outname);
};

template<class T>
dynDebug& dynDebug::operator<<(const T& t)
{
    if (status != ignore)
        *streamout << t;
    return *this;
}

#endif
