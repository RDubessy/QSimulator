/* Copyright (C) 2013 Romain Dubessy */
#ifndef GPE3D_H
#define GPE3D_H
#include "gpe.h"
#include "thermal.h"
/* class GPE3D {{{ */
/*!\brief This class implements a Gross Pitaevskii equation in a three
 * dimensional space.
 */
class GPE3D : public GPE {
    public:
        /*!\brief Constructor. */
        GPE3D(ConfigMap &config, Expression *H, Expression *pot);
        void spectrum(string &name, int m=0);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        double norm() const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
        void computePhase(std::complex<double> dt);
        void initialize(Expression *pot);
        double ekin();
        /*!\brief Initialize resources for the FFT. */
        virtual void initializeFFT();
    protected:
        double _dx;     //!<Grid step size along X
        double _dy;     //!<Grid step size along Y
        double _dz;     //!<Grid step size along Z
        int _nx;        //!<Number of grid points along X
        int _ny;        //!<Number of grid points along Y
        int _nz;        //!<Number of grid points along Z
        double _xmax;   //!<Maximum value of X on the grid
        double _ymax;   //!<Maximum value of Y on the grid
        double _zmax;   //!<Maximum value of Z on the grid
};
/* }}} */
/* class GPE3DROT {{{ */
/*!\brief This class implements a Gross Pitaevskii equation in a three 
 * dimensional space whithin a rotating frame.
 */
class GPE3DROT : public GPE3D {
    public:
        /*!\brief Constructor. */
        GPE3DROT(ConfigMap &config, Expression *H, Expression *pot);
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
/* class GPE3DThermal {{{ */
class GPE3DThermal : public GPE3D, public Thermal {
    public:
        /*!\brief Constructor. */
        GPE3DThermal(ConfigMap &config, Expression *H, Expression *pot, VarDef &params);
        void plot(int nmode, std::string &name);
        std::string measure();
        void setHeader(std::ofstream &file) const;
        state getHeader(std::ifstream &file);
        void save(std::string &name) const;
        void load(std::string &name);
        void doStep(std::complex<double> dt);
        void findGroundState(double dttest, double tol, double dttol, string &name, int verb=0);
        double thermalStep();
};
/* }}} */
#endif //GPE3D_H
/* gpe.h */
