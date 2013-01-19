/* Copyright (C) 2013 Romain Dubessy */
#ifndef GPE_H
#define GPE_H
#include <cvm.h>        //For cvm namespace
#include <fftw3.h>      //For fftw
#include <common.h>     //For ConfigMap
#include <expression.h> //For Expression, VarDef
const Constant zero(0);
const Constant one(1);
const double pi=acos(-1);
enum state {
    ok,error,converted
};
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
        GPE(Expression *H, Expression *pot);
        /*!\brief Computes the groundstate wave function of the system. */
        virtual void findGroundState(double dttest, double tol, double dttol, std::string &name, int verb=0);
        /*!\brief Computes the Bogolyubov spectrum of the system. */
        virtual void spectrum(std::string &name, int m=0);
        /*!\brief Normalize the groundstate wave function. */
        double normalize();
        /*!\brief Computes the norm of a real vector. */
        virtual double norm(cvm::rvector &psi) const=0;
        /*!\brief Computes the norm of a complex vector. */
        virtual double norm(cvm::cvector &psi) const=0;
        /*!\brief Implements the --plot option through gnuplot. */
        virtual void plot(int nmode, std::string &name) =0;
        /*!\brief Save the groundstate wave function to a file. */
        virtual void save(std::string &name) const;
        /*!\brief Load the groundstate wave function from a file. */
        virtual void load(std::string &name);
        /*!\brief Write the header of an output file. */
        virtual void setHeader(std::ofstream &file) const=0;
        /*!\brief Get the header from an input file. */
        virtual state getHeader(std::ifstream &file);
        /*!\brief Correct the hamiltonian for excited states. */
        virtual void correct(cvm::srmatrix &H, int m) {};
        /*!\brief Computes one step of evolution. */
        virtual void doStep(std::complex<double> dt);
        /*!\brief Computes the real time evolution of the system. */
        void evolve(double tstart, double dttest, double tend, std::string &name);
        /*!\brief Memory allocation method. */
        void allocate(int n);
        /*!\brief Compute the phases for Fourier methods. */
        virtual void computePhase(std::complex<double> dt) {};
        /*!\brief Initialize the system. */
        virtual void initialize(Expression *pot) =0;
        /*!\brief Evaluate the Hamiltonian. */
        virtual void update(double t);
        /*!\brief Perform a measurement on the system. */
        virtual std::string measure();
        /*!\brief Computes the potential energy. */
        virtual double epot();
        /*!\brief Computes the kinetic energy. */
        virtual double ekin() { return 0; };
        /*!\brief Imprint an integer winding number. */
        virtual void imprint(int l) {};
    protected:
        const char *_type;
        cvm::srbmatrix _H0; //!<Single body hamiltonian, stored as a tridiagonal real matrix.
        cvm::cvector _psi;  //!<Groundstate wave function, stored as a complex vector.
        double _gN;         //!<Interaction term.
        double _kterm;      //!<Kinetic term prefactor.
        double _vterm;      //!<Potential term prefactor.
        double _mu;         //!<Chemical potential (or groundstate energy if no interactions).
        std::complex<double> *_psip;    //!<Fourier space wave function representation.
        std::complex<double> *_phase;   //!<Phase for Fourier transform methods.
        double *_vpot;      //!<Real space potential array.
        fftw_plan _planFFT; //!<Resource for forward FFT.
        fftw_plan _planIFFT;//!<Resource for backward FFT.
        Expression *_H;     //!<Symbolic representation of the Hamiltonian.
    private:
        Expression *_pot;   //!<Symbolic representation of the potential.
};
/* }}} */
#endif //GPE_H
/* gpe.h */
