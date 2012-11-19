#ifndef MOVABLE_H
#define MOVABLE_H

class abstract_param {
    public: 
        virtual void generate_tweak(double theta_bar) = 0;
        virtual void restore_tweak() = 0;
};

class movable {
    public:
        double propose_move();
        void accept_move();
        void reject_move();
        void init_stats();
        virtual double get_score() = 0;
        virtual ~movable();
        movable();

    protected:
        int nparams;
        abstract_param **params;

    private:
        int index;
        long *success;
        long *moves;
        double *theta_bars;
        double energy;
        double prev_eng;
};


#endif
