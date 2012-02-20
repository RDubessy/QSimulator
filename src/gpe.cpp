#include <common.h>
#include <expression.h>
#include <fstream>
#include "gpe.h"
Constant zero(0);
Constant one(1);
double pi=acos(-1);
/* Quick Sort {{{ */
void swap(std::complex<double> &a, std::complex<double> &b) {
    std::complex<double> tmp=a;
    a=b;
    b=tmp;
};
int partition(std::complex<double> *v, int first, int last, int pivot,cvm::scmatrix &evecs) {
    swap(v[pivot],v[last]);
    evecs.swap_cols(pivot+1,last+1);
    int j=first;
    for(int i=first;i<last;i++) {
        if(v[i].real()<v[last].real()) {
            swap(v[i],v[j]);
            evecs.swap_cols(i+1,j+1);
            j++;
        }
    }
    swap(v[j],v[last]);
    evecs.swap_cols(j+1,last+1);
    return j;
};
void quickSort(std::complex<double> *v, int first, int last,cvm::scmatrix &evecs) {
    if(first<last) {
        int pivot=first+(rand()%(last-first+1));
        pivot=partition(v,first,last,pivot,evecs);
        quickSort(v,first,pivot-1,evecs);
        quickSort(v,pivot+1,last,evecs);
    }
    return;
};
/* }}} */
/* class GPE implementation {{{ */
/* findGroundState method {{{ */
void GPE::findGroundState(double dttest, double tol) {
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
}
/* }}} */
/* spectrum method {{{ */
void GPE::spectrum() {
    cvm::srbmatrix H1(_H);
    cvm::srbmatrix H3(_H);
    int n=_psi.size();
    double *v1=new double[n];
    double *v3=new double[n];
    double *v=_psi.get();
    for(int i=0;i<n;i++) {
        double gterm=_gN*v[i]*v[i];
        v1[i]=gterm-_mu;
        v3[i]=3*gterm-_mu;
    }
    H3.diag(0)+=cvm::rvector(v3,n);
    H1.diag(0)+=cvm::rvector(v1,n);
    delete[] v1;
    delete[] v3;
    cvm::srbmatrix H=H3*H1;
    cvm::scmatrix evecs(n);
    cvm::cvector evals(n);
    evals=H.eig(evecs);
    quickSort(evals.get(),0,n-1,evecs);
    for(int i=1;i<=n;i++) {
        std::cout << i << ' ' << evals(i).real()
            << ' ' << evals(i).imag() << '\n';
    }
    std::cout.flush();
    cvm::scmatrix upv=cvm::scmatrix(H1)*evecs;
    for(int i=1;i<=n;i++) {
        cvm::cvector u(n),v(n);
        u=(upv(i)/evals(i)+evecs(i))/2;
        v=(upv(i)/evals(i)-evecs(i))/2;
        double normu=norm(u);
        double normv=norm(v);
        double delta=normu*normu-normv*normv;
        u/=sqrt(delta);
        v/=sqrt(delta);
    }
}
/* }}} */
/* normalize method {{{ */
double GPE::normalize() {
    double res=norm(_psi);
    _psi/=res;
    return res;
};
/* }}} */
/* save method {{{ */
void GPE::save(std::string &name) const {
    std::ofstream file(name.c_str(),std::ofstream::binary);
    if(file.is_open()) {
        setHeader(file);
        const double *v=_psi.get();
        int n=_psi.size();
        for(int i=0;i<n;i++) {
            file.write((const char*)&(v[i]),sizeof(double));
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file:" << name << std::endl;
    }
    return;
}
/* }}} */
/* load method {{{ */
void GPE::load(std::string &name) {
    std::ifstream file(name.c_str(),std::ifstream::binary);
    if(file.is_open()) {
        if(!getHeader(file))
            return;
        double *v=_psi.get();
        int n=_psi.size();
        for(int i=0;i<n;i++) {
            file.read((char*)&(v[i]),sizeof(double));
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file:" << name << std::endl;
    }
    return;
}
/* }}} */
/* }}} */
/* class PolarGPE implementation {{{ */        
/* Constructor {{{ */
PolarGPE::PolarGPE(ConfigMap &config, Expression *H,
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
}
/* }}} */
/* norm methods {{{ */
double PolarGPE::norm(cvm::rvector &psi) const {
    double res=0.;
    double *v=psi.get();
    int n=psi.size();
    double rmin=_rmin/_dr;
    for(int i=0;i<n;i++)
        res+=(rmin+i)*v[i]*v[i];
    return _dr*sqrt(2*pi*res);
}
double PolarGPE::norm(cvm::cvector &psi) const {
    double res=0.;
    std::complex<double> *v=psi.get();
    int n=psi.size();
    double rmin=_rmin/_dr;
    for(int i=0;i<n;i++)
        res+=(rmin+i)*std::norm(v[i]);
    return _dr*sqrt(2*pi*res);
}
/* }}} */
/* plot method {{{ */
void PolarGPE::plot() const {
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
}
/* }}} */
void PolarGPE::setHeader(std::ofstream &file) const {
    const char *type="Pol";
    file.write((const char*)type,3*sizeof(char));
    file.write((const char*)&_n,sizeof(int));
    file.write((const char*)&_rmin,sizeof(double));
    file.write((const char*)&_dr,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
bool PolarGPE::getHeader(std::ifstream &file) {
    char type[3];
    double mu,dr,rmin;
    int n;
    file.read((char*)type,3*sizeof(char));
    file.read((char*)&n,sizeof(int));
    file.read((char*)&rmin,sizeof(double));
    file.read((char*)&dr,sizeof(double));
    file.read((char*)&mu,sizeof(double));
    _mu=mu;
    if(type[0]!='P' || type[1]!='o' || type[2]!='l') {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return false;
    }
    if(n!=_n) {
        std::cerr << "[E] Incompatible size !" << std::endl;
        return false;
    }
    std::cerr << "[I] Loading wavefunction (mu : " << _mu << ", #grid : " << _n << ")" << std::endl;
    return true;
}
/* }}} */
/* gpe.cpp */
