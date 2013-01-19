/* Copyright (C) 2013 Romain Dubessy */
#ifndef GPE2D_H
#define GPE2D_H
#include "gpe.h"
#include "thermal.h"
/* class Polar1D {{{ */
/*!\brief This class implements a Gross-Pitaevskii equation in polar coordinates
 * (2D space) assuming rotational invariance, which allows to reduce the problem
 * to an effective 1D equation.
 *
 * In this case the laplacian is:
 * \f[
 * \Delta=\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}
 * -\frac{\ell^2}{r^2},
 * \f]
 * where \f$\ell\f$ is a quantum number associated to the rotational invariance,
 * namely:
 * \f[
 * \hat{L}_z\psi=\hbar\ell\psi.
 * \f]
 */
class Polar1D : public GPE {
    public:
        /*!\brief Constructor. */
        Polar1D(ConfigMap &config, Expression *H, Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
        void correct(cvm::srmatrix &H, int m);
        void doStep(std::complex<double> dt);
        void findGroundState(double dttest, double tol, double dttol, string &name, int verb=0);
        void initialize(Expression *pot);
    protected:
        double _rmin;   //!<Minimum allowed radius.
        double _dr;     //!<Grid step size.
        int _n;         //!<Number of grid points.
        int _l;         //!<Quantum number associated to the rotationnal invariance.
    private:
        double _rmax;   //!<Maximum allowed radius.
};
/* }}} */
/* class Polar1DThermal {{{ */
class Polar1DThermal : public Polar1D, public Thermal {
    public:
        /*!\brief Constructor. */
        Polar1DThermal(ConfigMap &config, Expression *H, Expression *pot, VarDef &params);
        void doStep(std::complex<double> dt);
        std::string measure();
        void plot(int nmode, std::string &name);
        void findGroundState(double dttest, double tol, double dttol, string &name, int verb=0);
        double thermalStep();
        void save(std::string &name) const;
        void load(std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
};
/* }}} */
/* class GPE2D {{{ */
/*!\brief This class implements a Gross Pitaevskii equation in a two
 * dimensional space.
 *
 * In this case the laplacian is:
 * \f[\Delta=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}.\f]
 */
class GPE2D : public GPE {
    public:
        /*!\brief Constructor. */
        GPE2D(ConfigMap &config, Expression *H, Expression *pot);
        void spectrum(string &name, int m=0);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
        void computePhase(std::complex<double> dt);
        void initialize(Expression *pot);
        double ekin();
        void imprint(int l);
        /*!\brief Initialize resources for the FFT. */
        virtual void initializeFFT();
    protected:
        double _dx;     //!<Grid step size along X
        double _dy;     //!<Grid step size along Y
        int _nx;        //!<Number of grid points along X
        int _ny;        //!<Number of grid points along Y
        double _xmax;   //!<Maximum value of X on the grid
        double _ymax;   //!<Maximum value of Y on the grid
};
/* }}} */
/* class GPE2DROT {{{ */
/*!\brief This class implements a Gross Pitaevskii equation in a two dimensional
 * space whithin a rotating frame.
 *
 * In this case the laplacian is:
 * \f[\Delta=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}.\f]
 */
class GPE2DROT : public GPE2D {
    public:
        /*!\brief Constructor. */
        GPE2DROT(ConfigMap &config, Expression *H, Expression *pot);
        void doStep(std::complex<double> dt);
        void computePhase(std::complex<double> dt);
        void initializeFFT();
        void update(double dt);
        std::string measure();
    private:
        double _oterm;          //!<Rotation term (angular rotation frequency).
        fftw_plan _planFFTxz;   //!<Resource for forward FFT.
        fftw_plan _planFFTy;    //!<Resource for forward FFT.
        fftw_plan _planIFFTx;   //!<Resource for backward FFT.
        fftw_plan _planIFFTyz;  //!<Resource for backward FFT.
        std::complex<double> *_phase2;//!<Correction to the kinetic energy contribution.
};
/* }}} */
/* class GPE2DThermal {{{ */
class GPE2DThermal : public GPE2D, public Thermal {
    public:
        /*!\brief Constructor. */
        GPE2DThermal(ConfigMap &config, Expression *H, Expression *pot, VarDef &params);
        void doStep(std::complex<double> dt);
        std::string measure();
        void plot(int nmode, std::string &name);
        void findGroundState(double dttest, double tol, double dttol, string &name, int verb=0);
        double thermalStep();
        void save(std::string &name) const;
        void load(std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
};
/* }}} */
#endif //GPE2D_H
/* gpe.h */
