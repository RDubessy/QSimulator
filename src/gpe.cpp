#include <common.h>
#include <expression.h>
#include <fstream>
#include <iomanip>
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
/* Constructor {{{ */
GPE::GPE(Expression *H, Expression *pot) {
    _H=H;
    _pot=pot;
    _vpot=0;
    update(0);
    _psip=0;
    _phase=0;
    _vpot=0;
}
/* }}} */
/* findGroundState method {{{ */
/*! This method uses an imaginary time propagation algorithm to compute the
 * groundstate of the system:
 * \f[
 * \left|\psi(t+dt)\right>=\exp{\left[-\frac{Ht}{\hbar}\right]}\left|\psi(t)\right>.
 * \f]
 * After each step the norm of the function is computed which allows to evaluate
 * its chemical potential:
 * \f[
 * \mu=-\frac{\log{\left[\left<\psi(t+dt)\right|
 * \left.\psi(t+dt)\right>\right]}}{dt}.
 * \f]
 * The algorithm terminates when the relative change in the chemical potential
 * is below a given value.
 *
 * Since a high accuracy is needed in the determination of the groundstate, it
 * is more efficient to first choose the evolution time step be monitoring the
 * norm and when the norm change is sufficiently small keep the same time step
 * to find the groundstate.
 *
 * The initial state \f$\left|\psi(0)\right>\f$ must be initialized before calling
 * this function which allows to apply this algorithm to a guess groundstate that
 * must be refined.
 *
 * \param dttest an initial (largest) time step.
 * \param tol the relative target accuracy for the chemical potential.
 * \param dttol the relative target accuracy on the fonction norm.
 * \param name a string containing the name of a log file.
 */
void GPE::findGroundState(double dttest, double tol, double dttol, string &name) {
    std::cerr << "[I] Find groundstate method..." << std::endl;
    int c=0;
    int p=1;
    double dt=dttest;
    double muOld=1e3;
    double eps=1.0;
    double ttol=0.01;
    std::ofstream file(name.c_str());
    bool logout=file.is_open();
    if(!logout && name.size()>0) {
        cerr << "[W] Cannot open logging file: " << name << std::endl;
    }
    _mu=0;
    //Initialize the phase for fft
    computePhase(std::complex<double>(-dt,0));
    while(ttol>tol) {
        std::cerr << "[I]\tpass #" << p << " objective: " << ttol;
        std::cerr.flush();
        while(eps>ttol) {
            doStep(std::complex<double>(-dt,0));
            double tmp=normalize();
            if(tmp<dttol) {        //Choose dt.
                dt/=2;
                //Update the phase
                computePhase(std::complex<double>(-dt,0));
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
        //Update the phase
        computePhase(std::complex<double>(-dt,0));
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
    A+=_H0;
    D-=_H0;
    double *psi2=new double[n];
    double *d=new double[n];
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++) {
        psi2[i]=_gN*std::norm(v[i]);
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
                file.write((const char*)&(E[i-1]),sizeof(std::complex<double>));
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
        const std::complex<double> *v=_psi.get();
        int n=_psi.size();
        for(int i=0;i<n;i++) {
            file.write((const char*)&(v[i].real()),sizeof(double));
        }
        for(int i=0;i<n;i++) {
            file.write((const char*)&(v[i].imag()),sizeof(double));
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
        switch(getHeader(file)) {
            case ok:
                {
                    std::complex<double> *v=_psi.get();
                    int n=_psi.size();
                    std::cerr << "[I] Loading wavefunction (mu : " << _mu << ", #grid : " << n << ")" << std::endl;
                    for(int i=0;i<n;i++) {
                        file.read((char*)&(v[i].real()),sizeof(double));
                    }
                    if(!file.eof()) {
                        for(int i=0;i<n;i++) {
                            file.read((char*)&(v[i].imag()),sizeof(double));
                        }
                    }
                }
            case error:
            case converted:
            default:
                file.close();
        }
    } else {
        std::cerr << "[E] Can't open file:" << name << std::endl;
    }
    return;
}
/* }}} */
/* doStep method {{{ */
void GPE::doStep(std::complex<double> dt) {
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++)
        v[i]*=std::exp(dt*(_vpot[i]+_gN*std::norm(v[i])-_mu));
    fftw_execute(_planFFT);
    for(int i=0;i<n;i++)
        _psip[i]*=_phase[i];
    fftw_execute(_planIFFT);
}
/* }}} */
/* evolve method {{{ */
void GPE::evolve(double tstart, double dttest, double tend, std::string &name) {
    double t=tstart;
    std::complex<double> dt=std::complex<double>(0,-dttest);
    int c=0;
    update(t);
    computePhase(dt);
    while(t<tend) {
        doStep(dt);
        if(c%100==0) {
            std::cout << t << ' ' << norm(_psi)
                << ' ' << measure()
                << '\n';
        }
        if(c%1000==0) {
            std::cout.flush();
            std::ostringstream file;
            file << name << "_t" << std::setfill('0') << std::setw(4) << t;
            std::string filename=file.str();
            save(filename);
        }
        t+=dttest;
        //update(t);
        //computePhase(dt);
        c++;
    }
}
/* }}} */
/* allocate method {{{ */
void GPE::allocate(int n) {
    _psip=new std::complex<double>[n];
    _phase=new std::complex<double>[n];
    _vpot=new double[n];
}
/* }}} */
/* update method {{{ */
void GPE::update(double t) {
    VarDef vars;
    vars["LZ"]=&zero;
    vars["t_"]=new Constant(t);
    //Kinetic term pre-factor
    vars["DELTA"]=&one;
    vars["VEXT"]=&zero;
    vars["RHO"]=&zero;
    _kterm=*((double*)(_H->evaluate(vars)));
    //Interaction term pre-factor
    vars["DELTA"]=&zero;
    vars["VEXT"]=&zero;
    vars["RHO"]=&one;
    _gN=*((double*)(_H->evaluate(vars)));
    //Potential term pre-factor
    vars["DELTA"]=&zero;
    vars["VEXT"]=&one;
    vars["RHO"]=&zero;
    _vterm=*((double*)(_H->evaluate(vars)));
}

/* }}} */
/* measure method {{{ */
std::string GPE::measure() {
    std::ostringstream mes;
    mes << epot() << ' ' << ekin();
    std::string res=mes.str();
    return res;
}
/* }}} */
/* epot method {{{ */
double GPE::epot() {
    double res=0.;
    double n2=0.;
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++) {
        double v2=std::norm(v[i]);
        res+=v2*(_vpot[i]+_gN*v2/2);
        n2+=v2;
    }
    return res/n2;
}
/* }}} */
/* }}} */
/* class Polar1D implementation {{{ */ 
/* Constructor {{{ */
Polar1D::Polar1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 1D Polar system" << std::endl;
    _n=getConfig(config,string("polar::n"),64);
    _rmin=getConfig(config,string("polar::rmin"),0.01);
    _rmax=getConfig(config,string("polar::rmax"),1.0);
    _dr=(_rmax-_rmin)/(_n-1);
    _l=getConfig(config,string("polar::l"),0);
    _psi.resize(_n);
    _H0.resize(_n);
    _H0.resize_lu(1,1);
    allocate(_n);
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
    if(nmodes>0 && spectrum.is_open() && (getHeader(spectrum)==ok)) {
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
        const std::complex<double> *psi=_psi.get();
        for(int i=0;i<_n;i++) {
            file << _rmin+i*_dr << ' ' << psi[i].real() << ' ' << psi[i].imag();
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
        << "plot \"/tmp/psi.txt\" using 1:($2*$2+$3*$3) title \"\"";
    for(int i=0;i<nmodes;i++) {
        std::cout << ",\"\" using 1:(2*$2*($" << 4*i+4 << "+$" << 4*i+6 << ")) title \"\"";
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
    file.write((char*)&_l,sizeof(int));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state Polar1D::getHeader(std::ifstream &file) {
    char type[3];
    double mu,dr,rmin;
    int n,l;
    file.read((char*)type,3*sizeof(char));
    if(type[0]!='P' || type[1]!='o' || type[2]!='l') {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return error;
    }
    file.read((char*)&n,sizeof(int));
    file.read((char*)&rmin,sizeof(double));
    file.read((char*)&dr,sizeof(double));
    file.read((char*)&l,sizeof(int));
    file.read((char*)&mu,sizeof(double));
    _mu=mu;
    if(n!=_n) {
        std::cerr << "[E] Incompatible size !" << std::endl;
        return error;
    }
    return ok;
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
/* doStep method {{{ */
void Polar1D::doStep(std::complex<double> dt) {
    int n=_psi.size();
    cvm::cvector psi(n);
    std::complex<double> *v=psi.get();
    std::complex <double> *_v=_psi.get();
    psi=_psi;
    _psi.real()=_H0*psi.real();
    _psi.imag()=_H0*psi.imag();
    //Computes the interaction term
    for(int i=0;i<n;i++)
        _v[i]+=std::norm(v[i])*v[i]*_gN;
    _psi*=dt;
    _psi+=psi;
}
/* }}} */
/* initialize method {{{ */
void Polar1D::initialize(Expression *pot) {
    //Diagonal part
    VarDef vars;
    vars["R"]=new Constant(0);
    double *v=new double[_n];
    double *psi=new double[_n];
    for(int i=0;i<_n;i++) {
        double r=_rmin+i*_dr;
        vars["R"]->set(&r);
        double vpot=*((double*)(pot->evaluate(vars)));
        vpot*=_vterm;
        _vpot[i]=vpot;
        v[i]=_kterm*(-2/(_dr*_dr)-_l*_l/(r*r))+vpot;
        psi[i]=vpot<1?sqrt(1-vpot):0;
    }
    //Special case r=0
    vars["R"]->set(&_rmin);
    v[0]=*((double*)(pot->evaluate(vars)));
    v[0]*=_vterm;
    v[0]+=_kterm*(-1/(_dr*_dr)-_l*_l/(_rmin*_rmin));
    //Special case r=inf
    double rmax=_rmin+(_n-1)*_dr;
    vars["R"]->set(&rmax);
    v[_n-1]=*((double*)(pot->evaluate(vars)));
    v[_n-1]*=_vterm;
    v[_n-1]+=_kterm*(-1/(_dr*_dr)-_l*_l/(rmax*rmax));
    //Assign H diagonal and initial state
    _H0.diag(0).assign(v);
    _psi=cvm::cvector(psi,_n);
    delete[] v;
    //Upper and Lower diagonals
    double *vu=new double[_n-1];
    double *vl=new double[_n-1];
    for(int i=0;i<_n-1;i++) {
        double r=_rmin+i*_dr;
        vu[i]=_kterm*(1./_dr+0.5/r)/_dr;
        vl[i]=_kterm*(1./_dr-0.5/(r+_dr))/_dr;
    }
    //r=0
    vu[0]=_kterm/(_dr*_dr);
    _H0.diag(1).assign(vu);
    _H0.diag(-1).assign(vl);
    delete[] vu;
    delete[] vl;
}
/* }}} */
/* }}} */
/* class GPE1D implementation {{{ */        
/* Constructor {{{ */
GPE1D::GPE1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 1D Cartesian system" << std::endl;
    _n=getConfig(config,string("x1D::n"),64);
    _xmax=getConfig(config,string("x1D::L"),0.01);
    _dx=_xmax/(_n-1);
    _xmax/=2;
    _psi.resize(_n);
    _H0.resize(_n);
    _H0.resize_lu(1,1);
    allocate(_n);
}
/* }}} */
/* norm methods {{{ */
double GPE1D::norm(cvm::rvector &psi) const {
    return sqrt(_dx)*psi.norm2();
}
double GPE1D::norm(cvm::cvector &psi) const {
    return sqrt(_dx)*psi.norm2();
}
/* }}} */
/* plot method {{{ */
void GPE1D::plot(int nmodes, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    std::ifstream spectrum(name.c_str());
    std::complex<double> *u=0;
    std::complex<double> *v=0;
    if(nmodes>0 && spectrum.is_open() && (getHeader(spectrum)==ok)) {
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
        const std::complex<double> *psi=_psi.get();
        for(int i=0;i<_n;i++) {
            file << i*_dx-_xmax << ' ' << psi[i].real() << ' ' << psi[i].imag();
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
        << "plot \"/tmp/psi.txt\" using 1:($2*$2+$3*$3) title \"\"";
    for(int i=0;i<nmodes;i++) {
        std::cout << ",\"\" using 1:(2*$2*($" << 2*i+4 << "+$" << 2*i+5 << ")) title \"\"";
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
state GPE1D::getHeader(std::ifstream &file) {
    char type[3];
    double mu,dx,xmax;
    int n;
    file.read((char*)type,3*sizeof(char));
    if(type[0]!='x' || type[1]!='1' || type[2]!='D') {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return error;
    }
    file.read((char*)&n,sizeof(int));
    file.read((char*)&xmax,sizeof(double));
    file.read((char*)&dx,sizeof(double));
    file.read((char*)&mu,sizeof(double));
    _mu=mu;
    if(n!=_n) {
        std::cerr << "[E] Incompatible size !" << std::endl;
        return error;
    }
    return ok;
}
/* }}} */
/* computePhase method {{{ */
void GPE1D::computePhase(std::complex<double> dt) {
    double scale=1./_n;
    for(int i=0;i<_n;i++) {
        double k2=2*(cos((2*pi*i)/_n)-1)/(_dx*_dx);
        _phase[i]=scale*exp(dt*_kterm*k2);
    }
}
/* }}} */
/* initialize method {{{ */
void GPE1D::initialize(Expression *pot) {
    //Diagonal part
    VarDef vars;
    vars["X"]=new Constant(0);
    double *v=new double[_n];
    double *psi=new double[_n];
    for(int i=0;i<_n;i++) {
        double x=i*_dx-_xmax;
        vars["X"]->set(&x);
        double vpot=*((double*)(pot->evaluate(vars)));
        vpot*=_vterm;
        _vpot[i]=vpot;
        v[i]=_kterm*(-2/(_dx*_dx))+vpot;
        psi[i]=vpot<1?sqrt(1-vpot):0;
    }
    _H0.diag(0).assign(v);
    delete[] v;
    //Upper and lower diagonal
    double *vu=new double[_n-1];
    double *vl=new double[_n-1];
    for(int i=0;i<_n-1;i++) {
        vu[i]=_kterm/(_dx*_dx);
        vl[i]=vu[i];
    }
    _H0.diag(1).assign(vu);
    _H0.diag(-1).assign(vl);
    delete[] vu;
    delete[] vl;
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi.get());
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    _planFFT=fftw_plan_dft_3d(1,1,_n,rspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    _planIFFT=fftw_plan_dft_3d(1,1,_n,pspace,rspace,FFTW_BACKWARD,FFTW_MEASURE);
    //Copy initial state into memory
    _psi=cvm::cvector(psi,_n);
    delete[] psi;
}
/* }}} */
/* }}} */
/* class GPE2D implementation {{{ */
/* Constructor {{{ */
GPE2D::GPE2D(ConfigMap &config, Expression *H, Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 2D Cartesian system" << std::endl;
    _nx=getConfig(config,string("x2D::nx"),64);
    _ny=getConfig(config,string("x2D::ny"),_nx);
    _xmax=getConfig(config,string("x2D::Lx"),1.);
    _ymax=getConfig(config,string("x2D::Ly"),_xmax);
    _dx=_xmax/(_nx-1);
    _dy=_ymax/(_ny-1);
    _xmax/=2;
    _ymax/=2;
    int n=_nx*_ny;
    _psi.resize(n);
    _H0.resize(n);
    allocate(n);
}
/* }}} */
/* spectrum method {{{ */
void GPE2D::spectrum(string &name, int m) {
    std::cerr << "[E] spectrum method NOT implemented for a 2D system !" << std::endl;
    return;
}
/* }}} */
/* norm methods {{{ */
double GPE2D::norm(cvm::rvector &psi) const {
    return sqrt(_dx*_dy)*psi.norm2();
}
double GPE2D::norm(cvm::cvector &psi) const {
    return sqrt(_dx*_dy)*psi.norm2();
}
/* }}} */
/* plot method {{{ */
void GPE2D::plot(int nmode, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    std::ifstream spectrum(name.c_str());
    if(file.is_open()) {
        const std::complex<double> *psi=_psi.get();
        for(int j=0;j<_ny;j++) {
            for(int i=0;i<_nx;i++) {
                file << i*_dx-_xmax << ' ' << j*_dy-_ymax << ' '
                    << psi[i+j*_nx].real() << ' ' << psi[i+j*_nx].imag()
                    << '\n';
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi.txt\" !" << std::endl;
    }
    std::cout << "set view map;unset surface;set pm3d;unset key;"
        << "set size square;set size 1,1;set origin 0,0;"
        << "set multiplot layout 1,2;"
        << "splot \"/tmp/psi.txt\" using 1:2:($3*$3+$4*$4);"
        << "splot \"/tmp/psi.txt\" using 1:2:(arg($3+{0,1}*$4)+pi);"
        << "unset multiplot;"
        << "pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* setHeader method {{{ */
void GPE2D::setHeader(std::ofstream &file) const {
    const char *type="x2D";
    file.write((const char*)type,3*sizeof(char));
    file.write((const char*)&_nx,sizeof(int));
    file.write((const char*)&_ny,sizeof(int));
    file.write((const char*)&_xmax,sizeof(double));
    file.write((const char*)&_ymax,sizeof(double));
    file.write((const char*)&_dx,sizeof(double));
    file.write((const char*)&_dy,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state GPE2D::getHeader(std::ifstream &file) {
    char type[3];
    file.read((char*)type,3*sizeof(char));
    if(type[0]=='x' && type[1]=='2' && type[2]=='D') {
        double mu,dx,xmax,dy,ymax;
        int nx,ny;
        file.read((char*)&nx,sizeof(int));
        file.read((char*)&ny,sizeof(int));
        file.read((char*)&xmax,sizeof(double));
        file.read((char*)&ymax,sizeof(double));
        file.read((char*)&dx,sizeof(double));
        file.read((char*)&dy,sizeof(double));
        file.read((char*)&mu,sizeof(double));
        _mu=mu;
        if(nx!=_nx || ny!=_ny) {
            std::cerr << "[E] Incompatible size !" << std::endl;
            return error;
        }
    } else if(type[0]=='P' || type[1]=='o' || type[2]=='l') { //Convert from Polar1D
        std::cerr << "[I] Converting from polar coordinates" << std::endl;
        int n,l;
        double rmin,dr,mu;
        file.read((char*)&n,sizeof(int));
        file.read((char*)&rmin,sizeof(double));
        file.read((char*)&dr,sizeof(double));
        file.read((char*)&l,sizeof(int));
        file.read((char*)&mu,sizeof(double));
        std::complex<double> *v=new std::complex<double>[n];
        for(int i=0;i<n;i++) {
            file.read((char*)&(v[i].real()),sizeof(double));
        }
        if(!file.eof()) {
            for(int i=0;i<n;i++) {
                file.read((char*)&(v[i].imag()),sizeof(double));
            }
        }
        _mu=mu;
        std::complex<double> *w=_psi.get();
        double xo=0.5*(_nx-1);
        double yo=0.5*(_ny-1);
        double rmax=rmin+(n-1)*dr;
        for(int j=0;j<_ny;j++) {
            double y=_dy*(j-yo);
            for(int i=0;i<_nx;i++) {
                double x=_dx*(i-xo);
                double r=sqrt(x*x+y*y);
                if(r<rmin || r>rmax)
                    w[i+j*_nx]=0;
                else {
                    double theta=atan(y/x);
                    w[i+j*_nx]=v[(int)((r-rmin)/dr)]*std::exp(std::complex<double>(0,l*theta));
                }
            }
        }
        normalize();
        return converted;
    } else {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return error;
    }
    return ok;
}
/* }}} */
/* computePhase method {{{ */
void GPE2D::computePhase(std::complex<double> dt) {
    double scale=1./(_nx*_ny);
    for(int j=0;j<_ny;j++) {
        double k2y=(cos((2*pi*j)/_ny)-1)/(_dy*_dy);
        for(int i=0;i<_nx;i++) {
            double k2x=(cos((2*pi*i)/_nx)-1)/(_dx*_dx);
            double k2=2*(k2x+k2y);
            _phase[i+j*_nx]=scale*exp(dt*_kterm*k2);
        }
    }
}
/* }}} */
/* initialize method {{{ */
void GPE2D::initialize(Expression *pot) {
    //Only a diagonal part
    VarDef vars;
    vars["X"]=new Constant(0);
    vars["Y"]=new Constant(0);
    int n=_nx*_ny;
    double *v=new double[n];
    double *psi=new double[n];
    for(int j=0;j<_ny;j++) {
        double y=j*_dy-_ymax;
        vars["Y"]->set(&y);
        for(int i=0;i<_nx;i++) {
            double x=i*_dx-_xmax;
            vars["X"]->set(&x);
            double vpot=*((double*)(pot->evaluate(vars)));
            vpot*=_vterm;
            _vpot[i+j*_nx]=vpot;
            v[i+j*_nx]=vpot;
            psi[i+j*_nx]=vpot<1?sqrt(1-vpot):0;
        }
    }
    //Assign H diagonal and initial state
    _H0.diag(0).assign(v);
    delete[] v;
    initializeFFT();
    _psi=cvm::cvector(psi,n);
    delete[] psi;
}
/* }}} */
/* initializeFFT method {{{ */
void GPE2D::initializeFFT() {
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi.get());
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    _planFFT=fftw_plan_dft_3d(1,_ny,_nx,rspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    _planIFFT=fftw_plan_dft_3d(1,_ny,_nx,pspace,rspace,FFTW_BACKWARD,FFTW_MEASURE);
}
/* }}} */
/* ekin method {{{ */
double GPE2D::ekin() {
    double res=0.;
    double n2=0.;
    for(int j=0;j<_ny;j++) {
        double k2y=(cos((2*pi*j)/_ny)-1)/(_dy*_dy);
        for(int i=0;i<_nx;i++) {
            double k2x=(cos((2*pi*i)/_nx)-1)/(_dx*_dx);
            double psi2=std::norm(_psip[i+j*_nx]);
            res+=2*_kterm*(k2x+k2y)*psi2;
            n2+=psi2;
        }
    }
    return res/n2;
}
/* }}} */
/* }}} */
/* class GPE2DROT implementation {{{ */
/* Constructor {{{ */
GPE2DROT::GPE2DROT(ConfigMap &config, Expression *H, Expression *pot) : GPE2D(config,H,pot) {
    _phase2=new std::complex<double>[_nx*_ny];
    update(0);
}
/* }}} */
/* initializeFFT method {{{ */
void GPE2DROT::initializeFFT() {
    std::cerr << "[I] FFT in Rotating Frame" << std::endl;
    //Initialize the plans for fast fourier transform
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi.get());
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    fftw_iodim dimx,dimy,dimz;
    dimz.n=1;
    dimz.is=_nx*_ny;
    dimz.os=_nx*_ny;
    dimy.n=_ny;
    dimy.is=_nx;
    dimy.os=_nx;
    dimx.n=_nx;
    dimx.is=1;
    dimx.os=1;
    fftw_iodim dims[2];
    fftw_iodim hdims[2];
    //Transform x-z
    dims[0]=dimz;
    dims[1]=dimx;
    hdims[0]=dimy;
    _planFFTxz=fftw_plan_guru_dft(2,dims,1,hdims,rspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    //Transform y
    dims[0]=dimy;
    hdims[0]=dimz;
    hdims[1]=dimx;
    _planFFTy=fftw_plan_guru_dft(1,dims,2,hdims,pspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    //Inverse transform x
    dims[0]=dimx;
    hdims[0]=dimz;
    hdims[1]=dimy;
    _planIFFTx=fftw_plan_guru_dft(1,dims,2,hdims,pspace,pspace,FFTW_BACKWARD,FFTW_MEASURE);
    //Inverse transform y-z
    dims[0]=dimz;
    dims[1]=dimy;
    hdims[0]=dimx;
    _planIFFTyz=fftw_plan_guru_dft(2,dims,1,hdims,pspace,rspace,FFTW_BACKWARD,FFTW_MEASURE);
}
/* }}} */
/* computePhase method {{{ */
void GPE2DROT::computePhase(std::complex<double> dt) {
    double scale=1./(_nx*_ny);
    double xo=0.5*(_nx-1);
    double yo=0.5*(_ny-1);
    for(int j=0;j<_ny;j++) {
        double y=j-yo;
        double k2y=2*(cos((2*pi*j)/_ny)-1)/(_dy*_dy);
        double ky=sin((2*pi*j)/_ny);
        for(int i=0;i<_nx;i++) {
            double x=i-xo;
            double k2x=2*(cos((2*pi*i)/_nx)-1)/(_dx*_dx);
            double kx=sin((2*pi*i)/_nx);
            _phase[i+j*_nx]=scale*exp(dt*(_kterm*k2x-_oterm*y*kx));
            _phase2[i+j*_nx]=exp(dt*(_kterm*k2y+_oterm*x*ky));
        }
    }
    return;
}
/* }}} */
/* doStep method {{{ */
void GPE2DROT::doStep(std::complex<double> dt) {
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++)
        v[i]*=std::exp(dt*(_vpot[i]-_mu+_gN*std::norm(v[i])));
    fftw_execute(_planFFTxz);
    for(int i=0;i<n;i++)
        _psip[i]*=_phase[i];
    fftw_execute(_planFFTy);
    fftw_execute(_planIFFTx);
    for(int i=0;i<n;i++)
        _psip[i]*=_phase2[i];
    fftw_execute(_planIFFTyz);
}
/* }}} */
/* update method {{{ */
void GPE2DROT::update(double t) {
    GPE::update(t);
    VarDef vars;
    vars["t_"]=new Constant(t);
    //Rotation term pre-factor
    vars["LZ"]=&one;
    vars["DELTA"]=&zero;
    vars["VEXT"]=&zero;
    vars["RHO"]=&zero;
    _oterm=*((double*)(_H->evaluate(vars)));
}
/* }}} */
/* measure method {{{ */
std::string GPE2DROT::measure() {
    //Compute <Lz> & rotation energy
    double lz1=0.;
    double lz2=0.;
    double xo=0.5*(_nx-1);
    double yo=0.5*(_ny-1);
    double n2=0.;
    //Term ykx
    fftw_execute(_planFFTxz);
    for(int j=0;j<_ny;j++) {
        double y=j-yo;
        for(int i=0;i<_nx;i++) {
            double kx=sin((2*pi*i)/_nx);
            double psi2=std::norm(_psip[i+j*_nx]);
            lz1-=y*kx*psi2;
            n2+=psi2;
        }
    }
    lz1/=n2;
    fftw_execute(_planFFTy);
    //Compute kinetic energy: term kxky
    std::ostringstream mes;
    mes << GPE::measure();
    //Term xky
    n2=0.;
    fftw_execute(_planIFFTx);
    for(int j=0;j<_ny;j++) {
        double ky=sin((2*pi*j)/_ny);
        for(int i=0;i<_nx;i++) {
            double x=i-xo;
            double psi2=std::norm(_psip[i+j*_nx]);
            lz2+=x*ky*psi2;
            n2+=psi2;
        }
    }
    lz2/=n2;
    double lz=lz1+lz2;
    mes << ' ' << _oterm << ' ' << lz;
    std::string res=mes.str();
    return res;
}
/* }}} */
/* }}} */
/* gpe.cpp */
