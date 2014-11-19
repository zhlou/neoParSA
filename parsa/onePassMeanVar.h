/* 
 * File:   onePassMeanVar.h
 * Author: zhlou
 * 
 * A numerically stable one pass mean and variance algorithm a la TAOCP
 *
 * Created on November 19, 2014, 2:20 PM
 */

#ifndef ONEPASSMEANVAR_H
#define	ONEPASSMEANVAR_H


class onePassMeanVar {
private:
    double mean;
    double m2;
    size_t n;
public:
    onePassMeanVar() : mean(0), m2(0), n(0) {}
    void update(double data)
    {
        ++ n;
        double delta = data - mean;
        mean += delta / n;
        m2 += delta * (data - mean);
    }
    double getMean() const {return mean;}
    double getVar(size_t ddof=1) const {return (n<2) ? 0 : m2/(n-ddof);}
    void reset(){mean = 0; m2 = 0; n = 0;}
};



#endif	/* ONEPASSMEANVAR_H */

