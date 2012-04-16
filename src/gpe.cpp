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
void GPE::findGroundState(double dttest, double tol, double dttol, string &name) {
    std::cerr << "[I] Find groundstate method..." << std::endl;
    int c=0;
    int p=1;
    int n=_psi.size();
    double dt=dttest;
    double muOld=1e3;
    double eps=1.0;
    cvm::rvector psi(n);
    double *v=psi.get();
    double *_v=_psi.get();
    double ttol=0.01;
    std::ofstream file(name.c_str());
    bool logout=file.is_open();
    if(!logout && name.size()>0) {
        cerr << "[W] Cannot open logging file: " << name << std::endl;
    }
    while(ttol>tol) {
        std::cerr << "[I]\tpass #" << p << " objective: " << ttol;
        std::cerr.flush();
        while(eps>ttol) {
            psi=_psi;
            _psi=_H*psi;
            //Computes the interaction term
            for(int i=0;i<n;i++)
                _v[i]+=v[i]*v[i]*v[i]*_gN;
            _psi*=(-dt);
            _psi+=psi;
            double tmp=normalize();
            if(tmp<dttol) {        //Choose dt.
                dt/=2;
            } else {
                double mu=-log(tmp)/dt;
                eps=std::abs((mu-muOld)/muOld);
                muOld=mu;
            }
            c++;
            if(logout && c%100==0)
                file << c << ' ' << tmp << ' ' << dt << ' ' << muOld << ' ' << eps << '\n';
        }
        std::cerr << ", done: " << eps << ", mu=" << muOld << std::endl;
        ttol/=10;
        p++;
        dt=dttest;
    }
    if(logout)
        file.close();
    _mu=muOld;
    std::cerr << "[I] After "<< c << " iterations, mu=" << _mu << " [" << eps << "]." << std::endl;
    return;
}
/* }}} */
/* spectrum method {{{ */
void GPE::spectrum(string &name, int m) {
    std::cerr << "[I] Compute spectrum method... m=" << m << std::endl;
    int n=_psi.size();
    int N=2*n;
    cvm::srmatrix H(N);
    cvm::srmatrix A(H,1,1,n);
    cvm::srmatrix B(H,1,n+1,n);
    cvm::srmatrix C(H,n+1,1,n);
    cvm::srmatrix D(H,n+1,n+1,n);
    A+=_H;
    D-=_H;
    double *psi2=new double[n];
    double *d=new double[n];
    double *v=_psi.get();
    for(int i=0;i<n;i++) {
        psi2[i]=_gN*v[i]*v[i];
        d[i]=2*psi2[i]-_mu;
    }
    A.diag(0)+=cvm::rvector(d,n);
    D.diag(0)-=cvm::rvector(d,n);
    B.diag(0)+=cvm::rvector(psi2,n);
    C.diag(0)-=cvm::rvector(psi2,n);
    delete[] psi2;
    delete[] d;
    correct(H,m);
    cvm::cvector evals(N);
    cvm::scmatrix evecs(N);
    evals=H.eig(evecs);
    quickSort(evals.get(),0,N-1,evecs);
    char buffer[256];
    sprintf(buffer,"_spectrum_m%d",m);
    string s=name;
    s+=buffer;
    std::ofstream file(s.c_str(),std::ofstream::binary);
    double sumv2=0.;
    bool *pfamily=new bool[N];
    if(file.is_open()) {
        std::complex<double> *E=evals.get();
        setHeader(file);
        for(int i=1;i<=N;i++) {
            cvm::cvector col(evecs(i));
            cvm::cvector u(n),v(n);
            u=cvm::cvector(&(col.get()[0]),n);
            v=cvm::cvector(&(col.get()[n]),n);
            double normu=norm(u);
            double normv=norm(v);
            double delta=normu*normu-normv*normv;
            if(delta>0) {
                u/=sqrt(delta);
                v/=sqrt(delta);
                sumv2+=normv*normv/delta;
                std::complex<double> *U=u.get();
                std::complex<double> *V=v.get();
                file.write((const char*)&(E[i]),sizeof(std::complex<double>));
                for(int j=0;j<n;j++)
                    file.write((const char*)&(U[j]),sizeof(std::complex<double>));
                for(int j=0;j<n;j++)
                    file.write((const char*)&(V[j]),sizeof(std::complex<double>));
                pfamily[i-1]=true;
            } else
                pfamily[i-1]=false;
        }
    } else {
        std::cerr << "[E] Can't open file: " << s << std::endl;
    }
    if(name.size()>0) {
        std::ofstream logfile(name.c_str(),std::ofstream::app);
        bool logout=logfile.is_open();
        if(!logout) {
            cerr << "[W] Cannot open logging file: " << name << std::endl;
        } else {
            logfile << m << ' ' << sumv2;
            for(int i=1;i<=N;i++) {
                if(pfamily[i-1])
                    logfile << ' ' << evals(i).real()
                        << ' ' << evals(i).imag();
            }
            logfile << std::endl;
            logfile.close();
        }
    }
    return;
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
        std::cerr << "[I] Loading wavefunction (mu : " << _mu << ", #grid : " << n << ")" << std::endl;
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
/* class Polar1D implementation {{{ */        
/* Constructor {{{ */
Polar1D::Polar1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE() {
    _n=getConfig(config,string("polar::n"),64);
    _rmin=getConfig(config,string("polar::rmin"),0.01);
    _rmax=getConfig(config,string("polar::rmax"),1.0);
    _dr=(_rmax-_rmin)/(_n-1);
    _l=getConfig(config,string("polar::l"),0);
    _psi.resize(_n);
    _H.resize(_n);
    _H.resize_lu(1,1);
    VarDef vars;
    //Kinetic term pre-factor
    vars["DELTA"]=&one;
    vars["VEXT"]=&zero;
    vars["RHO"]=&zero;
    _kterm=*((double*)(H->evaluate(vars)));
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
    double r0=(_rmin+_rmax)/2;
    double R=sqrt(2*5);
    for(int i=0;i<_n;i++) {
        double r=_rmin+i*_dr;
        vars["R"]->set(&r);
        v[i]=*((double*)(diag->evaluate(vars)));
        //psi[i]=std::exp(-5*(*((double*)(pot->evaluate(vars))))/(_rmax-_rmin));
        psi[i]=std::abs(r-r0)<R?sqrt(5/_gN*(1-pow((r-r0)/R,2))):0;
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
double Polar1D::norm(cvm::rvector &psi) const {
    double res=0.;
    double *v=psi.get();
    int n=psi.size();
    double rmin=_rmin/_dr;
    for(int i=0;i<n;i++)
        res+=(rmin+i)*v[i]*v[i];
    return _dr*sqrt(2*pi*res);
}
double Polar1D::norm(cvm::cvector &psi) const {
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
void Polar1D::plot(int nmodes, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    std::ifstream spectrum(name.c_str());
    std::complex<double> *u=0;
    std::complex<double> *v=0;
    if(nmodes>0 && spectrum.is_open() && getHeader(spectrum)) {
        std::cerr << "[I] Found matching spectrum file" << std::endl;
        u=new std::complex<double>[_n*_n];
        v=new std::complex<double>[_n*_n];
        for(int i=0;i<_n;i++) {
            std::complex<double> e;
            spectrum.read((char *)&e,sizeof(std::complex<double>));
            for(int j=0;j<_n;j++)
                spectrum.read((char *)&u[i*_n+j],sizeof(std::complex<double>));
            for(int j=0;j<_n;j++)
                spectrum.read((char *)&v[i*_n+j],sizeof(std::complex<double>));
        }
        spectrum.close();
    } else
        nmodes=0;
    if(file.is_open()) {
        const double *psi=_psi.get();
        for(int i=0;i<_n;i++) {
            file << _rmin+i*_dr << ' ' << psi[i];
            for(int j=0;j<nmodes;j++) {
                file << ' ' << u[j*_n+i].real() << ' ' << u[j*_n+i].imag() << ' '
                    << ' ' << v[j*_n+i].real() << ' ' << v[j*_n+i].imag();
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi.txt\" !" << std::endl;
    }
    std::cout << "set style data lines;"
        << "set xlabel \"r\";set ylabel \"Density\";"
        << "plot \"/tmp/psi.txt\" using 1:($2*$2) title \"\"";
    for(int i=0;i<nmodes;i++) {
        std::cout << ",\"\" using 1:(2*$2*($" << 4*i+3 << "+$" << 4*i+5 << ")) title \"\"";
    }
    std::cout << ";pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* setHeader method {{{ */
void Polar1D::setHeader(std::ofstream &file) const {
    const char *type="Pol";
    file.write((const char*)type,3*sizeof(char));
    file.write((const char*)&_n,sizeof(int));
    file.write((const char*)&_rmin,sizeof(double));
    file.write((const char*)&_dr,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
bool Polar1D::getHeader(std::ifstream &file) {
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
    return true;
}
/* }}} */
/* correct method {{{ */
void Polar1D::correct(cvm::srmatrix &H, int m) {
    if(m==0) return;
    double *v=new double[2*_n];
    double cor=-1.*_kterm;
    for(int i=0;i<_n;i++) {
        double r=_rmin+i*_dr;
        double invr2=cor/(r*r);
        v[i]=m*(m+2*_l)*invr2;
        v[i+_n]=m*(2*_l-m)*invr2;
    }
    H.diag(0)+=cvm::rvector(v,2*_n);
    delete[] v;
}
/* }}} */
/* }}} */
/* class GPE1D implementation {{{ */        
/* Constructor {{{ */
GPE1D::GPE1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE() {
    _n=getConfig(config,string("x1D::n"),64);
    _xmax=getConfig(config,string("x1D::L"),0.01);
    _dx=_xmax/(_n-1);
    _xmax/=2;
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
    vars["DX"]=new Constant(_dx);
    vars["VEXT"]=pot;
    vars["DELTA"]=parseString("-2/DX^2");
    Expression *diag=H->simplify(vars);
    vars["VEXT"]=&zero;
    vars["DELTA"]=parseString("1/DX^2");
    Expression *diagu=H->simplify(vars);
    vars["DELTA"]=parseString("1/DX^2");
    Expression *diagl=H->simplify(vars);
    vars["X"]=new Constant(0);
    double *v=new double[_n];
    double *psi=new double[_n];
    for(int i=0;i<_n;i++) {
        double x=i*_dx-_xmax;
        vars["X"]->set(&x);
        v[i]=*((double*)(diag->evaluate(vars)));
        psi[i]=std::exp(-2.5*(*((double*)(pot->evaluate(vars))))/(_xmax));
    }
    _H.diag(0).assign(v);
    _psi.assign(psi);
    delete[] v;
    double *vu=new double[_n-1];
    double *vl=new double[_n-1];
    for(int i=0;i<_n-1;i++) {
        vu[i]=*((double*)(diagu->evaluate(vars)));
        vl[i]=*((double*)(diagl->evaluate(vars)));
    }
    _H.diag(1).assign(vu);
    _H.diag(-1).assign(vl);
    _H(1,1)=-vu[0];
    _H(_n,_n)=-vu[0];
    delete[] vu;
    delete[] vl;
}
/* }}} */
/* norm methods {{{ */
double GPE1D::norm(cvm::rvector &psi) const {
    double res=0.;
    double *v=psi.get();
    int n=psi.size();
    for(int i=0;i<n;i++)
        res+=v[i]*v[i];
    return sqrt(_dx*res);
}
double GPE1D::norm(cvm::cvector &psi) const {
    double res=0.;
    std::complex<double> *v=psi.get();
    int n=psi.size();
    for(int i=0;i<n;i++)
        res+=std::norm(v[i]);
    return sqrt(_dx*res);
}
/* }}} */
/* plot method {{{ */
void GPE1D::plot(int nmodes, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    std::ifstream spectrum(name.c_str());
    std::complex<double> *u=0;
    std::complex<double> *v=0;
    if(nmodes>0 && spectrum.is_open() && getHeader(spectrum)) {
        std::cerr << "[I] Found matching spectrum file" << std::endl;
        u=new std::complex<double>[_n*_n];
        v=new std::complex<double>[_n*_n];
        for(int i=0;i<_n;i++) {
            std::complex<double> e;
            spectrum.read((char *)&e,sizeof(std::complex<double>));
            for(int j=0;j<_n;j++)
                spectrum.read((char *)&u[i*_n+j],sizeof(std::complex<double>));
            for(int j=0;j<_n;j++)
                spectrum.read((char *)&v[i*_n+j],sizeof(std::complex<double>));
        }
        spectrum.close();
    } else
        nmodes=0;
    if(file.is_open()) {
        const double *psi=_psi.get();
        for(int i=0;i<_n;i++) {
            file << i*_dx-_xmax << ' ' << psi[i];
            for(int j=0;j<nmodes;j++) {
                file << ' ' << u[j*_n+i].real()
                    << ' ' << v[j*_n+i].real();
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi.txt\" !" << std::endl;
    }
    std::cout << "set style data lines;"
        << "set xlabel \"x\";set ylabel \"Density\";"
        << "plot \"/tmp/psi.txt\" using 1:($2*$2) title \"\"";
    for(int i=0;i<nmodes;i++) {
        std::cout << ",\"\" using 1:(2*$2*($" << 2*i+3 << "+$" << 2*i+4 << ")) title \"\"";
        //std::cout << ",\"\" using 1:($" << 2*i+4 << ") title \"\"";
    }
    std::cout << ";pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* setHeader method {{{ */
void GPE1D::setHeader(std::ofstream &file) const {
    const char *type="x1D";
    file.write((const char*)type,3*sizeof(char));
    file.write((const char*)&_n,sizeof(int));
    file.write((const char*)&_xmax,sizeof(double));
    file.write((const char*)&_dx,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
bool GPE1D::getHeader(std::ifstream &file) {
    char type[3];
    double mu,dx,xmax;
    int n;
    file.read((char*)type,3*sizeof(char));
    file.read((char*)&n,sizeof(int));
    file.read((char*)&xmax,sizeof(double));
    file.read((char*)&dx,sizeof(double));
    file.read((char*)&mu,sizeof(double));
    _mu=mu;
    if(type[0]!='x' || type[1]!='1' || type[2]!='D') {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return false;
    }
    if(n!=_n) {
        std::cerr << "[E] Incompatible size !" << std::endl;
        return false;
    }
    return true;
}
/* }}} */
/* }}} */
/* gpe.cpp */
