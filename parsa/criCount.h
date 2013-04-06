/*
 * criCount.h
 *
 *  Created on: Apr 5, 2013
 *      Author: zhlou
 */

#ifndef CRICOUNT_H_
#define CRICOUNT_H_

class criCountP;
class criCount {
public:
    class Param
    {
        double freeze_crit;
        int cnt_crit;
        Param(xmlNode *root);
    };
    criCount(const Param &param);
    ~criCount();
    void updateStep(bool, const aState &){}
    bool frozen(const aState &state);
    friend class criCountP;
private:
    const double freeze_crit;
    double old_energy;
    const int cnt_crit;
    int freeze_cnt;
};

#endif /* CRICOUNT_H_ */
