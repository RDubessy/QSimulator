/* Copyright (C) 2013 Romain Dubessy */
#ifndef THERMAL_H
#define THERMAL_H
double polylog3half(double z);
/*!\brief This class describe a semi-classical Bose gas.
 */
class Thermal {
    public:
        /*!\brief Constructor. */
        Thermal(ConfigMap &config, VarDef &params);
        /*!\brief Destructor. */
        ~Thermal();
        /*!\brief Computes the equilibrium thermal atoms density.*/
        virtual double thermalStep() =0;
    protected:
        double *_n0;    //!<Thermal atoms density, stored as a real array.
        double _Ntot;   //!<Total number of atoms.
        double _beta;   //!<Inverse of the temperature.
        double _lambda; //!<Thermal de Broglie wavelength.
        double _Nbec;   //!<Number of condensed atoms.
        double _Ntherm; //!<Number of thermal atoms.
};
#endif //THERMAL_H
/* thermal.h */
