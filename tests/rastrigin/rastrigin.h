#ifndef RASTRIGIN_H
#define RASTRIGIN_H
class variable : public abstract_param {
    public:
        static const VAR_MAX = 5.12;
        static const VAR_MIN = -5.12;
        double x;
        unsigned int *seed;
        void generate_tweek(double theta_bar);
        void restore_tweak();
}

class rastrigin : public movable {
    public:
        rastrigin (int dimension);
        ~rastrigin();
        double value() const;
        void init_vars(unsigned int *seed);
    private:
        variable *vars;

}

#endif
