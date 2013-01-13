/* Copyright (C) 2013 Romain Dubessy */
#ifndef THERMAL_H
#define THERMAL_H
class Thermal {
    public:
        /*!\brief Constructor. */
        Thermal(ConfigMap &config, VarDef &params);
        /*!\brief Destructor. */
        ~Thermal();
    protected:
        double *_n0;
        double _Ntot;
        double _beta;
        double _lambda;
        double _Nbec;
        double _Ntherm;
};
#endif //THERMAL_H
/* thermal.h */
