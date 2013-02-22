/* Copyright (C) 2013 Romain Dubessy */
#ifndef GPE1D_H
#define GPE1D_H
#include "gpe.h"
/* class GPE1D {{{ */
/*!\brief This class implements a Gross-Pitaevskii equation in a one dimensional
 * space.
 *
 * In this case the laplacian is:
 * \f[\Delta=\frac{\partial^2}{\partial x^2}.\f]
 */
class GPE1D : public GPE {
    public:
        /*!\brief Constructor. */
        GPE1D(ConfigMap &config, Expression *H, Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        double norm() const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
        void computePhase(std::complex<double> dt);
        void initialize(Expression *pot);
        double ekin();
    private:
        double _xmax;   //!<Half box size.
        double _dx;     //!<Grid step size.
        int _n;         //!<Number of grid points.
};
/* }}} */
#endif //GPE1D_H
/* gpe.h */
