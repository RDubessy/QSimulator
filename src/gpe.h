#ifndef GPE_H
#define GPE_H
#include <cvm.h>
#include <fftw3.h>
/* class GPE {{{ */
/*!\brief This class describe a Gross-Pitaevskii equation.
 *
 * This equation may be written in physical units:
 * \f[
 * \imath\hbar\partial_t\psi(r,t)=\left(-\frac{\hbar^2}{2m}\Delta+V(r)
 * +gN\left|\psi(r,t)\right|^2\right)\psi(r,t),
 * \f]
 * where \f$\hbar\f$ is the reduced Planck constant, \f$m\f$ is the particle
 * mass, \f$\Delta\f$ is the Laplacian, \f$V(r)\f$ is the external trapping
 * potential, \f$g\f$ is the interaction coupling constant and \f$N\f$ the atom
 * number.
 * This equation assumes the normalization condition:
 * \f[
 * \int d^3r \left|\psi(r,t)\right|^2=1.
 * \f]
 *
 * The main purpose is to compute the groundstate of the equation and the linear
 * spectrum of excitation (through the Bogolyubov de Gennes equations).
 * These tasks are done in to public methods:
 * - findGroundState(),
 * - spectrum().
 *
 * This class is virtual and therefore cannot be instancied, because the
 * Laplacian is not specified.
 */
class GPE {
    public:
        /*!\brief Constructor. */
        GPE(Expression *H);
        /*!\brief Computes the groundstate wave function of the system. */
        void findGroundState(double dttest, double tol, double dttol, string &name);
        /*!\brief Computes the Bogolyubov spectrum of the system. */
        virtual void spectrum(string &name, int m=0);
        /*!\brief Normalize the groundstate wave function. */
        double normalize();
        /*!\brief Computes the norm of a real vector. */
        virtual double norm(cvm::rvector &psi) const=0;
        /*!\brief Computes the norm of a complex vector. */
        virtual double norm(cvm::cvector &psi) const=0;
        /*!\brief Implements the --plot option through gnuplot. */
        virtual void plot(int nmode, std::string &name) =0;
        /*!\brief Save the groundstate wave function to a file. */
        void save(std::string &name) const;
        /*!\brief Load the groundstate wave function from a file. */
        void load(std::string &name);
        /*!\brief Write the header of an output file. */
        virtual void setHeader(std::ofstream &file) const=0;
        /*!\brief Get the header from an input file. */
        virtual bool getHeader(std::ifstream &file) =0;
        /*!\brief Correct the hamiltonian for excited states. */
        virtual void correct(cvm::srmatrix &H, int m) {};
        virtual void doStep(std::complex<double> dt);
        void evolve(double tstart, double dttest, double tend);
        void allocate(int n);
        virtual void computePhase(std::complex<double> dt) {};
        virtual void initialize(Expression *pot) =0;
    protected:
        cvm::srbmatrix _H;  //!<Single body hamiltonian, stored as a tridiagonal real matrix.
        cvm::cvector _psi;  //!<Groundstate wave function, stored as a complex vector.
        double _gN;         //!<Interaction term.
        double _kterm;      //!<Kinetic term prefactor.
        double _vterm;      //!<Potential term prefactor.
        double _mu;         //!<Chemical potential (or groundstate energy if no interactions).
        std::complex<double> *_psip;
        std::complex<double> *_phase;
        double *_vpot;
        fftw_plan _planFFT; //!<Resource for forward FFT.
        fftw_plan _planIFFT;//!<Resource for backward FFT.
};
/* }}} */
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
        Polar1D(ConfigMap &config, Expression *H,
                Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
        void correct(cvm::srmatrix &H, int m);
        void doStep(std::complex<double> dt);
        void initialize(Expression *pot);
    private:
        double _rmin;   //!<Minimum allowed radius.
        double _rmax;   //!<Maximum allowed radius.
        double _dr;     //!<Grid step size.
        int _l;         //!<Quantum number associated to the rotationnal invariance.
        int _n;         //!<Number of grid points.
};
/* }}} */
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
        GPE1D(ConfigMap &config, Expression *H,
                Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
        void computePhase(std::complex<double> dt);
        void initialize(Expression *pot);
    private:
        double _xmax;   //!<Half box size.
        double _dx;     //!<Grid step size.
        int _n;         //!<Number of grid points.
};
/* }}} */
/* class GPE2D {{{ */
class GPE2D : public GPE {
    public:
        GPE2D(ConfigMap &config, Expression *H,
                Expression *pot);
        void spectrum(string &name, int m=0);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode, std::string &name);
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
        void computePhase(std::complex<double> dt);
        void initialize(Expression *pot);
        virtual void initializeFFT();
    protected:
        double _dx;
        double _dy;
        int _nx;
        int _ny;
    private:
        double _xmax;
        double _ymax;
};
/* }}} */
/* class GPE2DROT {{{ */
class GPE2DROT : public GPE2D {
    public:
        GPE2DROT(ConfigMap &config, Expression *H,
                Expression *pot);
        void doStep(std::complex<double> dt);
        void computePhase(std::complex<double> dt);
        void initializeFFT();
    private:
        double _oterm;          //!<Rotation term.
        fftw_plan _planFFTxz;   //!<Resource for forward FFT.
        fftw_plan _planFFTy;    //!<Resource for forward FFT.
        fftw_plan _planIFFTx;   //!<Resource for backward FFT.
        fftw_plan _planIFFTyz;  //!<Resource for backward FFT.
        std::complex<double> *_phase2;//!<Correction to the kinetic energy contribution.
};
/* }}} */
#endif //GPE_H
/* gpe.h */
