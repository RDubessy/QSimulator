/* This file is a part of QSimulator. {{{
 * Copyright (C) 2013 Romain Dubessy
 *
 * QSimulator is free software: you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *
 * QSimulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QSimulator.  If not, see <http://www.gnu.org/licenses/>.
 *
 * }}} */
#include "gpe2D.h"
#include <fstream>      //For ifstream, ofstream
/* class Polar1D implementation {{{ */ 
/* Constructor {{{ */
Polar1D::Polar1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 1D Polar system" << std::endl;
    _type="pol";
    _n=getConfig(config,string("grid::nr"),64);
    _rmin=getConfig(config,string("grid::rmin"),0.01);
    _rmax=getConfig(config,string("grid::rmax"),1.0);
    _dr=(_rmax-_rmin)/(_n-1);
    _l=getConfig(config,string("grid::l"),0);
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
    file.write((const char*)_type,3*sizeof(char));
    file.write((const char*)&_n,sizeof(int));
    file.write((const char*)&_rmin,sizeof(double));
    file.write((const char*)&_dr,sizeof(double));
    file.write((char*)&_l,sizeof(int));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state Polar1D::getHeader(std::ifstream &file) {
    if(GPE::getHeader(file)==error)
        return error;
    int n;
    file.read((char*)&n,sizeof(int));
    file.read((char*)&_rmin,sizeof(double));
    file.read((char*)&_dr,sizeof(double));
    file.read((char*)&_l,sizeof(int));
    file.read((char*)&_mu,sizeof(double));
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
        _v[i]+=_gN*std::norm(v[i])*v[i];
    _psi*=dt;
    _psi+=psi;
}
/* }}} */
/* findGroundState method {{{ */
void Polar1D::findGroundState(double dttest, double tol, double dttol, string &name, int verb) {
    std::cerr << "[I] Find groundstate method..." << std::endl;
    int c=0;
    int p=1;
    double dt=dttest;
    double eps=1.0;
    double ttol=0.01;
    std::ofstream file(name.c_str());
    bool logout=file.is_open();
    if(!logout && name.size()>0) {
        cerr << "[W] Cannot open logging file: " << name << std::endl;
    }
    _mu=0;
    double muOld=0;
    //Initialize the phase for fft
    computePhase(std::complex<double>(-dt,0));
    while(ttol>tol) {
        if(verb>0)
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
                _mu=-std::log(tmp)/dt;
                eps=std::abs((muOld-_mu)/_mu);
                muOld=_mu;
            }
            c++;
            if(logout && c%100==0)
                file << c << ' ' << tmp << ' ' << dt << ' ' << _mu << ' ' << eps << '\n';
        }
        if(verb>0)
            std::cerr << ", done: " << eps << ", mu=" << _mu << std::endl;
        ttol/=10;
        p++;
        dt=dttest;
        //Update the phase
        computePhase(std::complex<double>(-dt,0));
    }
    if(logout)
        file.close();
    std::cerr << "[I] After "<< c << " iterations, mu=" << _mu << " [" << eps << "]." << std::endl;
    return;
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
/* class Polar1DThermal implementation {{{ */
/* Constructor {{{ */
Polar1DThermal::Polar1DThermal(ConfigMap &config, Expression *H, Expression *pot, VarDef &params) : Polar1D(config,H,pot), Thermal(config,params) {
    int n=_psi.size();
    _n0=new double[n];
    memset(_n0,0,n*sizeof(double));
    _type="P1T";
}
/* }}} */
/* doStep method {{{ */
void Polar1DThermal::doStep(std::complex<double> dt) {
    int n=_psi.size();
    cvm::cvector psi(n);
    std::complex<double> *v=psi.get();
    std::complex <double> *_v=_psi.get();
    psi=_psi;
    _psi.real()=_H0*psi.real();
    _psi.imag()=_H0*psi.imag();
    //Computes the interaction term
    for(int i=0;i<n;i++)
        _v[i]+=_gN*(_Nbec*std::norm(v[i])+2*_n0[i])*v[i];
    _psi*=dt;
    _psi+=psi;
}
/* }}} */
/* measure method {{{ */
std::string Polar1DThermal::measure() {
    std::ostringstream mes;
    mes.precision(15);
    mes << _mu << ' ' << epot() << ' ' << ekin() << ' ' << _Nbec << ' ' << _Ntherm;
    std::string res=mes.str();
    return res;
}
/* }}} */
/* plot method {{{ */
void Polar1DThermal::plot(int nmodes, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    if(file.is_open()) {
        const std::complex<double> *psi=_psi.get();
        for(int i=0;i<_n;i++) {
            file << _rmin+i*_dr << ' ' << psi[i].real() << ' ' << psi[i].imag() << ' ' << _n0[i] << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi.txt\" !" << std::endl;
    }
    std::cout << "set style data lines;"
        << "set xlabel \"r\";set ylabel \"Density\";"
        << "plot \"/tmp/psi.txt\" using 1:(($2*$2+$3*$3)*" << _Nbec << ") title \"\""
        << ",\"\" using 1:4 title \"\"";
    std::cout << ";pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* thermalStep method {{{ */
double Polar1DThermal::thermalStep() {
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++) {
        double z=std::exp(-_beta*(_vpot[i]+2*_gN*(_Nbec*std::norm(v[i])+_n0[i])-_mu));
        if(z<1)
            _n0[i]=-std::log(1.0-z)/(_lambda*_lambda);
    }
    double rmin=_rmin/_dr;
    double nth=rmin*_n0[0]+(rmin+n-1)*_n0[n-1];
    for(int i=1;i<n-1;i++) {
        _n0[i]=(_n0[i]+_n0[i-1])/2;
        nth+=(rmin+i)*_n0[i];
    }
    nth*=2*pi*_dr*_dr;
    return nth;
}
/* }}} */
/* findGroundState method {{{ */
void Polar1DThermal::findGroundState(double dttest, double tol, double dttol, string &name, int verb) {
    double nold,eps;
    int n=_psi.size();
    do {
        nold=_Ntherm;
        double ntmp,rel,mutmp,murel;
        do {
            ntmp=_Ntherm;
            mutmp=_mu;
            Polar1D::findGroundState(dttest,tol,dttol,name,verb);
            _Ntherm=thermalStep();
            rel=std::abs((ntmp-_Ntherm)/_Ntherm);
            murel=std::abs((mutmp-_mu)/_mu);
        } while(rel>1e-3 || murel>1e-7);
        double ntot=_Nbec+_Ntherm;
        double scale=_Ntot/ntot;
        _Nbec*=scale;
        _Ntherm*=scale;
        for(int i=0;i<n;i++)
            _n0[i]*=scale;
        std::cerr << "[I] Thermal step: [" << _Nbec << "/" << _Ntherm << "] " << std::endl;
        eps=std::abs(nold-_Ntherm);
    } while(eps>1);
}
/* }}} */
/* setHeader method {{{ */
void Polar1DThermal::setHeader(std::ofstream &file) const {
    Polar1D::setHeader(file);
    file.write((const char*)&_Nbec,sizeof(double));
    file.write((const char*)&_Ntherm,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state Polar1DThermal::getHeader(std::ifstream &file) {
    if(Polar1D::getHeader(file)==error)
        return error;
    file.read((char*)&_Nbec,sizeof(double));
    file.read((char*)&_Ntherm,sizeof(double));
    return ok;
}
/* }}} */
/* save method {{{ */
void Polar1DThermal::save(std::string &name) const {
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
        for(int i=0;i<n;i++) {
            file.write((const char*)&(_n0[i]),sizeof(double));
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file:" << name << std::endl;
    }
    return;
}
/* }}} */
/* load method {{{ */
void Polar1DThermal::load(std::string &name) {
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
                    for(int i=0;i<n;i++) {
                        file.read((char*)&(v[i].imag()),sizeof(double));
                    }
                    std::cerr << "[I] Loading thermal part [" << _Nbec << "/" << _Ntherm << "]" << std::endl;
                    for(int i=0;i<n;i++) {
                        file.read((char*)&(_n0[i]),sizeof(double));
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
/* }}} */
/* class GPE2D implementation {{{ */
/* Constructor {{{ */
GPE2D::GPE2D(ConfigMap &config, Expression *H, Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 2D Cartesian system" << std::endl;
    _type="x2D";
    _nx=getConfig(config,string("grid::nx"),64);
    _ny=getConfig(config,string("grid::ny"),_nx);
    _xmax=getConfig(config,string("grid::Lx"),1.);
    _ymax=getConfig(config,string("grid::Ly"),_xmax);
    _dx=_xmax/(_nx-1);
    _dy=_ymax/(_ny-1);
    _xmax/=2;
    _ymax/=2;
    int n=_nx*_ny;
    _psi.resize(n);
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
    file.write((const char*)_type,3*sizeof(char));
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
    if(type[0]==_type[0] && type[1]==_type[1] && type[2]==_type[2]) {
        int nx,ny;
        file.read((char*)&nx,sizeof(int));
        file.read((char*)&ny,sizeof(int));
        file.read((char*)&_xmax,sizeof(double));
        file.read((char*)&_ymax,sizeof(double));
        file.read((char*)&_dx,sizeof(double));
        file.read((char*)&_dy,sizeof(double));
        file.read((char*)&_mu,sizeof(double));
        if(nx!=_nx || ny!=_ny) {
            std::cerr << "[E] Incompatible size !" << std::endl;
            return error;
        }
    } else if(type[0]=='P' || type[1]=='o' || type[2]=='l') { //Convert from Polar1D
        std::cerr << "[I] Converting from polar coordinates" << std::endl;
        int n,l;
        double rmin,dr;
        file.read((char*)&n,sizeof(int));
        file.read((char*)&rmin,sizeof(double));
        file.read((char*)&dr,sizeof(double));
        file.read((char*)&l,sizeof(int));
        file.read((char*)&_mu,sizeof(double));
        std::complex<double> *v=new std::complex<double>[n];
        for(int i=0;i<n;i++) {
            file.read((char*)&(v[i].real()),sizeof(double));
        }
        if(!file.eof()) {
            for(int i=0;i<n;i++) {
                file.read((char*)&(v[i].imag()),sizeof(double));
            }
        }
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
        double ky=2*pi/(_dy*_ny)*(j<_ny/2?j:j-_ny);
        double k2y=ky*ky;
        for(int i=0;i<_nx;i++) {
            double kx=2*pi/(_dx*_nx)*(i<_nx/2?i:i-_nx);
            double k2x=kx*kx;
            double k2=-(k2x+k2y);
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
            psi[i+j*_nx]=vpot<1?sqrt(1-vpot):0;
        }
    }
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
/* imprint method {{{ */
void GPE2D::imprint(int l) {
    double xo=0.5*(_nx-1);
    double yo=0.5*(_ny-1);
    std::complex<double> *psi=_psi.get();
    for(int j=0;j<_ny;j++) {
        double y=j-yo;
        for(int i=0;i<_nx;i++) {
            double x=i-xo;
            double phase=-atan2(y,-x)+2*pi;
            psi[i+j*_nx]*=std::exp(std::complex<double>(0,l*phase));
        }
    }
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
        double ky=2*pi/(_dy*_ny)*(j<_ny/2?j:j-_ny);
        double k2y=-ky*ky;
        for(int i=0;i<_nx;i++) {
            double x=i-xo;
            double kx=2*pi/(_dx*_nx)*(i<_nx/2?i:i-_nx);
            double k2x=-kx*kx;
            _phase[i+j*_nx]=scale*std::exp(dt*(_kterm*k2x-_oterm*y*kx));
            _phase2[i+j*_nx]=std::exp(dt*(_kterm*k2y+_oterm*x*ky));
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
    vars["LZ"]=(Constant*)&one;
    vars["DELTA"]=(Constant*)&zero;
    vars["VEXT"]=(Constant*)&zero;
    vars["RHO"]=(Constant*)&zero;
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
/* class GPE2DThermal implementation {{{ */
/* Constructor {{{ */
GPE2DThermal::GPE2DThermal(ConfigMap &config, Expression *H, Expression *pot,VarDef &params) : GPE2D(config,H,pot), Thermal(config,params) {
    int n=_psi.size();
    _n0=new double[n];
    memset(_n0,0,n*sizeof(double));
    _type="T2D";
}
/* }}} */
/* doStep method {{{ */
void GPE2DThermal::doStep(std::complex<double> dt) {
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++)
        v[i]*=std::exp(dt*(_vpot[i]+_gN*(_Nbec*std::norm(v[i])+2*_n0[i])-_mu));
    fftw_execute(_planFFT);
    for(int i=0;i<n;i++)
        _psip[i]*=_phase[i];
    fftw_execute(_planIFFT);
}
/* }}} */
/* thermalStep method {{{ */
double GPE2DThermal::thermalStep() {
    int n=_psi.size();
    std::complex<double> *v=_psi.get();
    for(int i=0;i<n;i++) {
        double z=std::exp(-_beta*(_vpot[i]+2*_gN*(_Nbec*std::norm(v[i])+_n0[i])-_mu));
        if(z<1)
            _n0[i]=-std::log(1.0-z)/(_lambda*_lambda);
    }
    double nth=0;
    for(int i=0;i<n;i++) {
        nth+=_n0[i];
    }
    nth*=_dx*_dy;
    return nth;
}
/* }}} */
/* findGroundState method {{{ */
void GPE2DThermal::findGroundState(double dttest, double tol, double dttol, string &name, int verb) {
    double nold,eps;
    int n=_psi.size();
    do {
        nold=_Ntherm;
        double ntmp,rel,mutmp,murel;
        do {
            ntmp=_Ntherm;
            mutmp=_mu;
            GPE::findGroundState(dttest,tol,dttol,name,verb);
            _Ntherm=thermalStep();
            rel=std::abs((ntmp-_Ntherm)/_Ntherm);
            murel=std::abs((mutmp-_mu)/_mu);
        } while(rel>1e-3 || murel>1e-7);
        double ntot=_Nbec+_Ntherm;
        double scale=_Ntot/ntot;
        _Nbec*=scale;
        _Ntherm*=scale;
        for(int i=0;i<n;i++)
            _n0[i]*=scale;
        std::cerr << "[I] Thermal step: [" << _Nbec << "/" << _Ntherm << "] " << std::endl;
        eps=std::abs(nold-_Ntherm);
    } while(eps>1);
}
/* }}} */
/* measure method {{{ */
std::string GPE2DThermal::measure() {
    std::ostringstream mes;
    mes.precision(15);
    mes << _mu << ' ' << epot() << ' ' << ekin() << ' ' << _Nbec << ' ' << _Ntherm;
    std::string res=mes.str();
    return res;
}
/* }}} */
/* setHeader method {{{ */
void GPE2DThermal::setHeader(std::ofstream &file) const {
    GPE2D::setHeader(file);
    file.write((const char*)&_Nbec,sizeof(double));
    file.write((const char*)&_Ntherm,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state GPE2DThermal::getHeader(std::ifstream &file) {
    if(GPE2D::getHeader(file)==error)
        return error;
    file.read((char*)&_Nbec,sizeof(double));
    file.read((char*)&_Ntherm,sizeof(double));
    return ok;
}
/* }}} */
/* plot method {{{ */
void GPE2DThermal::plot(int nmode, std::string &name) {
    std::ofstream file("/tmp/psi.txt");
    if(file.is_open()) {
        const std::complex<double> *psi=_psi.get();
        for(int j=0;j<_ny;j++) {
            for(int i=0;i<_nx;i++) {
                file << i*_dx-_xmax << ' ' << j*_dy-_ymax << ' '
                    << psi[i+j*_nx].real() << ' ' << psi[i+j*_nx].imag()
                    << ' ' << _n0[i+j*_nx]
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
        << "splot \"/tmp/psi.txt\" using 1:2:5;"
        << "unset multiplot;"
        << "pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* save method {{{ */
void GPE2DThermal::save(std::string &name) const {
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
        for(int i=0;i<n;i++) {
            file.write((const char*)&(_n0[i]),sizeof(double));
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file:" << name << std::endl;
    }
    return;
}
/* }}} */
/* load method {{{ */
void GPE2DThermal::load(std::string &name) {
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
                    for(int i=0;i<n;i++) {
                        file.read((char*)&(v[i].imag()),sizeof(double));
                    }
                    std::cerr << "[I] Loading thermal part [" << _Nbec << "/" << _Ntherm << "]" << std::endl;
                    for(int i=0;i<n;i++) {
                        file.read((char*)&(_n0[i]),sizeof(double));
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
/* }}} */
/* gpe.cpp */
