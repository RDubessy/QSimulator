#ifndef GPE_H
#define GPE_H
#include <cvm.h>
/* class GPE {{{ */
class GPE {
    public:
        GPE() {};
        void findGroundState(double dttest, double tol);
        void spectrum(void);
        double normalize();
        virtual double norm(cvm::rvector &psi) const=0;
        virtual double norm(cvm::cvector &psi) const=0;
        virtual void plot() const=0;
        void save(std::string &name) const;
        void load(std::string &name);
        virtual void setHeader(std::ofstream &file) const=0;
        virtual bool getHeader(std::ifstream &file) =0;
    protected:
        cvm::srbmatrix _H;
        cvm::rvector _psi;
        double _gN;
        double _mu;
};
/* }}} */
/* class PolarGPE {{{ */
class PolarGPE : public GPE {
    public:
        PolarGPE(ConfigMap &config, Expression *H,
                Expression *pot);
        double norm(cvm::rvector &psi) const;
        double norm(cvm::cvector &psi) const;
        void plot() const;
        void setHeader(std::ofstream &file) const;
        bool getHeader(std::ifstream &file);
    private:
        double _rmin;
        double _rmax;
        double _dr;
        int _l;
        int _n;
};
/* }}} */
#endif //GPE_H
/* gpe.h */
