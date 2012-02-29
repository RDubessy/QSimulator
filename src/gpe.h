#ifndef GPE_H
#define GPE_H
#include <cvm.h>
/* class GPE {{{ */
class GPE {
    public:
        GPE() {};
        void findGroundState(double dttest, double tol);
        void spectrum(int m=0);
        double normalize();
        virtual double norm(cvm::rvector &psi) const=0;
        virtual double norm(cvm::cvector &psi) const=0;
        virtual void plot(int nmode) =0;
        void save(std::string &name) const;
        void load(std::string &name);
        virtual void setHeader(std::ofstream &file) const=0;
        virtual bool getHeader(std::ifstream &file) =0;
        virtual void correct(cvm::srbmatrix &H, int m) {};
    protected:
        cvm::srbmatrix _H;  //!<Single body hamiltonian, stored as a tridiagonal real matrix.
        cvm::rvector _psi;  //!<Groundstate wave function, stored as a real vector.
        double _gN;         //!<Interaction term.
        double _mu;         //!<Chemical potential (or groundstate energy if no interactions).
};
/* }}} */
/* class PolarGPE {{{ */
class PolarGPE : public GPE {
    public:
        PolarGPE(ConfigMap &config, Expression *H,
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
#endif //GPE_H
/* gpe.h */
