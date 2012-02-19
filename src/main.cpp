#include <iostream>     //For cerr/cout/endl...
#include <fstream>      //For ofstream/ifstream...
#include <common.h>     //For ConfigMap.
#include <expression.h> //For Expression, VarDef.
#include <cvm.h>        //For cvm matrices.
Constant zero(0);
Constant one(1);
double pi=acos(-1);
class GPE {
    public:
        GPE() {};
        /* findGroundState method {{{ */
        void findGroundState(double dttest, double tol) {
            int c=0;
            int n=_psi.size();
            double dt=dttest;
            double muOld=1e3;
            double eps=1.0;
            cvm::rvector psi(n);
            double *v=psi.get();
            double *_v=_psi.get();
            while(eps>tol) {
                psi=_psi;
                _psi=_H*psi;
                //Computes the interaction term
                for(int i=0;i<n;i++)
                    _v[i]+=v[i]*v[i]*v[i]*_gN;
                _psi*=(-dt);
                _psi+=psi;
                double tmp=normalize();
                if(tmp<0.999) {
                    dt/=2;
                } else {
                    double mu=-log(tmp)/dt;
                    eps=std::abs((mu-muOld)/muOld);
                    muOld=mu;
                }
                c++;
                if(c%100==0)
                    std::cout << c << ' ' << tmp << ' ' << dt << ' ' << muOld << ' ' << eps << '\n';
            }
            std::cout.flush();
            _mu=muOld;
            std::cerr << c << ' ' << norm(_psi) << ' ' << _mu << ' ' << eps << std::endl;
            return;
        };
        /* }}} */
        /* normalize method {{{ */
        double normalize() {
            double res=norm(_psi);
            _psi/=res;
            return res;
        };
        /* }}} */
        virtual double norm(cvm::rvector &psi) const=0;
        virtual double norm(cvm::cvector &psi) const=0;
        virtual void save(std::string &name) const=0;
        virtual void load(std::string &name)=0;
        virtual void plot() const=0;
    protected:
        cvm::srbmatrix _H;
        cvm::rvector _psi;
        double _gN;
        double _mu;
};
class PolarGPE : public GPE {
    public:
        /* Constructor {{{ */
        PolarGPE(ConfigMap &config, Expression *H,
                Expression *pot) : GPE() {
            _n=getConfig(config,string("polar::n"),64);
            _rmin=getConfig(config,string("polar::rmin"),0.01);
            _rmax=getConfig(config,string("polar::rmax"),1.0);
            _dr=(_rmax-_rmin)/_n;
            _l=0;
            _psi.resize(_n);
            _H.resize(_n);
            _H.resize_lu(1,1);
            VarDef vars;
            //Interaction term pre-factor
            vars["DELTA"]=&zero;
            vars["VEXT"]=&zero;
            vars["RHO"]=&one;
            _gN=*((double*)(H->evaluate(vars)));
            //Diagonal part
            vars["RHO"]=&zero;
            vars["DR"]=new Constant(_dr);
            vars["L"]=new Constant(_l);
            vars["VEXT"]=pot;
            vars["DELTA"]=parseString("-2/DR^2-(L/R)^2");
            Expression *diag=H->simplify(vars);
            vars["VEXT"]=&zero;
            vars["DELTA"]=parseString("1/DR^2+0.5/(R*DR)");
            Expression *diagu=H->simplify(vars);
            vars["DELTA"]=parseString("1/DR^2-0.5/((R+DR)*DR)");
            Expression *diagl=H->simplify(vars);
            vars["R"]=new Constant(0);
            double *v=new double[_n];
            double *psi=new double[_n];
            for(int i=0;i<_n;i++) {
                double r=_rmin+i*_dr;
                vars["R"]->set(&r);
                v[i]=*((double*)(diag->evaluate(vars)));
                psi[i]=std::exp(-5*(*((double*)(pot->evaluate(vars))))/(_rmax-_rmin));
            }
            _H.diag(0).assign(v);
            _psi.assign(psi);
            delete[] v;
            double *vu=new double[_n-1];
            double *vl=new double[_n-1];
            for(int i=0;i<_n-1;i++) {
                double r=_rmin+i*_dr;
                vars["R"]->set(&r);
                vu[i]=*((double*)(diagu->evaluate(vars)));
                vl[i]=*((double*)(diagl->evaluate(vars)));
            }
            _H.diag(1).assign(vu);
            _H.diag(-1).assign(vl);
            delete[] vu;
            delete[] vl;
        };
        /* }}} */
        /* norm methods {{{ */
        double norm(cvm::rvector &psi) const {
            double res=0.;
            double *v=psi.get();
            int n=psi.size();
            double rmin=_rmin/_dr;
            for(int i=0;i<n;i++)
                res+=(rmin+i)*v[i]*v[i];
            return _dr*sqrt(2*pi*res);
        };
        double norm(cvm::cvector &psi) const {
            double res=0.;
            std::complex<double> *v=psi.get();
            int n=psi.size();
            double rmin=_rmin/_dr;
            for(int i=0;i<n;i++)
                res+=(rmin+i)*std::norm(v[i]);
            return _dr*sqrt(2*pi*res);
        };
        /* }}} */
        /* save method {{{ */
        void save(std::string &name) const {
            std::ofstream file(name.c_str(),std::ofstream::binary);
            if(file.is_open()) {
                const char *type="Pol";
                file.write((const char*)type,3*sizeof(char));
                file.write((const char*)&_n,sizeof(int));
                file.write((const char*)&_rmin,sizeof(double));
                file.write((const char*)&_dr,sizeof(double));
                file.write((const char*)&_mu,sizeof(double));
                const double *v=_psi.get();
                for(int i=1;i<=_n;i++) {
                    file.write((const char*)&(v[i]),sizeof(double));
                }
                file.close();
            } else {
                std::cerr << "[E] Can't open file:" << name << std::endl;
            }
            return;
        };
        /* }}} */
        /* load method {{{ */
        void load(std::string &name) {
            std::ifstream file(name.c_str(),std::ifstream::binary);
            if(file.is_open()) {
                char type[3];
                double mu,dr,rmin;
                int n;
                file.read((char*)type,3*sizeof(char));
                file.read((char*)&n,sizeof(int));
                file.read((char*)&rmin,sizeof(double));
                file.read((char*)&dr,sizeof(double));
                file.read((char*)&mu,sizeof(double));
                if(type[0]!='P' || type[1]!='o' || type[2]!='l') {
                    std::cerr << "[E] Incompatible type !" << std::endl;
                    return;
                }
                if(n!=_n) {
                    std::cerr << "[E] Incompatible size !" << std::endl;
                    return;
                }
                std::cerr << "[I] Loading wavefunction (mu : " << mu << ", #grid : " << _n << ")" << std::endl;
                _mu=mu;
                double *v=_psi.get();
                for(int i=1;i<=n;i++) {
                    file.read((char*)&(v[i]),sizeof(double));
                }
                file.close();
            } else {
                std::cerr << "[E] Can't open file:" << name << std::endl;
            }
            return;
        };
        /* }}} */
        /* plot method {{{ */
        void plot() const {
            std::ofstream file("/tmp/psi.txt");
            if(file.is_open()) {
                const double *v=_psi.get();
                for(int i=0;i<_n;i++) {
                    file << _rmin+i*_dr << ' ' << v[i] << '\n';
                }
                file.close();
            } else {
                std::cerr << "[E] Can't open file \"/tmp/psi.txt\" !" << std::endl;
            }
            std::cout << "set xlabel \"r\";set ylabel \"Density\";"
                << "plot \"/tmp/psi.txt\" using 1:($2*$2) title \"\";pause mouse\n";
            std::cout.flush();
            return;
        };
        /* }}} */
    private:
        double _rmin;
        double _rmax;
        double _dr;
        int _l;
        int _n;
};
/* mainFunction method {{{ */
int mainFunction(ConfigMap &config) {
    if(config["general::equation"].size()==0) {
        cerr << "[E] No problem defined!" << endl;
        return -1;
    }
    cerr << "[I] Initializing problem...\n";
    Expression *eqn=parseString(config["general::equation"]);
    Expression *pot=parseString(config["general::potential"]);
    VarDef params;
    cerr << "[I] Getting parameters:\n";
    for(ConfigMap::iterator it=config.begin();it!=config.end();it++) {
        if((it->first).find("parameters::")!=string::npos) {
            string name=(it->first).substr(12);
            params[name]=parseString(it->second);
            cerr << "[I]\t" << name << "=" << it->second << "\n";
        }
    }
    cerr.flush();
    eqn=eqn->simplify(params);
    pot=pot->simplify(params);
    cerr << "[I] Equation: " << eqn << endl;
    cerr << "[I] Potential: " << pot << endl;
    GPE *gpe=0;
    if(config["general::type"]=="Polar") {
        gpe=new PolarGPE(config,eqn,pot);
    }
    if(gpe==0) {
        std::cerr << "[E] Unknown problem type!" << std::endl;
        return -1;
    }
    if(config["in"].size()>0) {
        gpe->load(config["in"]);
    }
    string out="psi1.dat";
    if(config["out"].size()>0) {
        out=config["out"];
    }
    if(config["groundstate"].size()>0) {
        double dt=getConfig(config,string("general::dt"),1e-3);
        double tol=getConfig(config,string("general::tol"),1e-8);
        gpe->findGroundState(dt,tol);
        gpe->save(out);
    }
    if(config["plot"].size()>0) {
        gpe->plot();
    }
    return 0;
};
/* }}} */
/* usage method {{{ */
void usage(const char *name) {
    cerr << name << endl;
};
/* }}} */
/* main method {{{ */
int main(int argc, char *argv[]) {
    ConfigMap config;
    if(!parseOptions(argc,argv,config)) {        //Parse cmd line options
        cerr << "==> Try '" << argv[0] << " --usage'" << endl;
        return -1;
    }
    if(config["usage"].size()>0||config["help"].size()>0) {
        usage(argv[0]);
        return -1;
    }
    if(!parseConfig(config)) {
        cerr << "==> Try 'man " << argv[0] << "'" << endl;
        return -1;
    }
    return mainFunction(config);
}
/* }}} */
/* main.cpp */
