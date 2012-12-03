#ifndef MOVABLE_H
#define MOVABLE_H
#include <libxml/parser.h>
class abstract_param {
    public: 
        virtual void generate_tweak(double theta_bar) = 0;
        virtual void restore_tweak() = 0;
        virtual ~abstract_param(); // need this on obj w/ virtual func.
};

class movable {
    public:
        double propose_move();
        void accept_move();
        void reject_move();
        void init_stats();
        void set_theta(int id, double theta);
        virtual double get_score() = 0;
        virtual ~movable();
        movable(int np);
        movable(xmlNode *root, xmlChar *secname);

    protected:
        int nparams;
        abstract_param **params;
        xmlNode *docroot;
        xmlNode *section;

    private:
        void init(int np);
        int index;
        long *success;
        long *moves;
        double *theta_bars;
        double energy;
        double prev_eng;
};


#endif
