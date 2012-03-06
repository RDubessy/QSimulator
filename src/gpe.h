#ifndef GPE_H
#define GPE_H
#include <cvm.h>
/* class GPE {{{ */
/*!\brief This class describe a Gross-Pitaevskii equation.
 *
 * This equation may be written in physical units:
 * \f[
 * \imath\partial_t\psi(r,t)=\left(-\frac{\hbar^2}{2m}\Delta+V(r)
 * +gN\left|\psi(r,t)\right|^2\right)\psi(r,t),
 * \f]
 * where \f$\hbar\f$ is the reduced Planck constant, \f$m\f$ is the particle
 * mass, \f$\Delta\f$ is the Laplacian, \f$V(r)\f$ is the external trapping
 * potential, \f$g\f$ is the interaction coupling constant and \f$N\f$ the atom
 * number.
 *
 * The main purpose is to compute the groundstate of the equation and the linear
 * spectrum of excitation (through the Bogolyubov de Gennes equations).
 *
 * This class is virtual and therefore cannot be instancied, because the
 * Laplacian is not specified.
 */
class GPE {
    public:
        /*!\brief Constructor. */
        GPE() {};
        /*!\brief Computes the groundstate wave function of the system. */
        void findGroundState(double dttest, double tol);
        /*!\brief Computes the Bogolyubov spectrum of the system. */
        void spectrum(int m=0);
        /*!\brief Normalize the groundstate wave function. */
        double normalize();
        /*!\brief Computes the norm of a real vector. */
        virtual double norm(cvm::rvector &psi) const=0;
        /*!\brief Computes the norm of a complex vector. */
        virtual double norm(cvm::cvector &psi) const=0;
        /*!\brief Implements the --plot option through gnuplot. */
        virtual void plot(int nmode) =0;
        /*!\brief Save the groundstate wave function to a file. */
        void save(std::string &name) const;
        /*!\brief Load the groundstate wave function from a file. */
        void load(std::string &name);
        /*!\brief Write the header of an output file. */
        virtual void setHeader(std::ofstream &file) const=0;
        /*!\brief Get the header from an input file. */
        virtual bool getHeader(std::ifstream &file) =0;
        /*!\brief Correct the hamiltonian for excited states. */
        virtual void correct(cvm::srbmatrix &H, int m) {};
    protected:
        cvm::srbmatrix _H;  //!<Single body hamiltonian, stored as a tridiagonal real matrix.
        cvm::rvector _psi;  //!<Groundstate wave function, stored as a real vector.
        double _gN;         //!<Interaction term.
        double _mu;         //!<Chemical potential (or groundstate energy if no interactions).
};
/* }}} */
/* class PolarGPE {{{ */
/*!\brief This class implements a Gross-Pitaevskii equation in polar coordinates
 * assuming rotational invariance.
 *
 * In this case the laplacian is:
 * \f[
 * \Delta=\frac{\partial^2}{\partial r^2}+\frac{1}{r}\frac{\partial}{\partial r}
 * -\frac{l^2}{r^2},
 * \f]
 * where \f$l\f$ is a quantum number associated to the rotational invariance,
 * namely:
 * \f[
 * \hat{L}_z\psi=\hbar l\psi.
 * \f]
 */
class Polar1D : public GPE {
    public:
        /*!\brief Constructor. */
        Polar1D(ConfigMap &config, Expression *H,
                Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot(int nmode);
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
        void correct(cvm::srbmatrix &H, int m);
    private:
        double _rmin;   //!<Minimum allowed radius.
        double _rmax;   //!<Maximum allowed radius.
        double _dr;     //!<Grid step size.
        double _kterm;  //!<Kinetic term prefactor.
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
        void plot(int nmode);
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
    private:
        double _xmax;   //!<Half box size.
        double _dx;     //!<Grid step size.
        int _n;         //!<Number of grid points.
};
/* }}} */
#endif //GPE_H
/* gpe.h */
