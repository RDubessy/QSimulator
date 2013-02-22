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
#include "gpe1D.h"
#include <fstream>      //For ifstream, ofstream
/* class GPE1D implementation {{{ */        
/* Constructor {{{ */
GPE1D::GPE1D(ConfigMap &config, Expression *H,
        Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 1D Cartesian system" << std::endl;
    _type="x1D";
    _n=getConfig(config,string("grid::nx"),64);
    _xmax=getConfig(config,string("grid::Lx"),0.01);
    _dx=_xmax/(_n-1);
    _xmax/=2;
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
double GPE1D::norm() const {
    double res=0;
    for(int i=0;i<_size;i++)
        res+=std::norm(_psi[i]);
    res*=_dx;
    return sqrt(res);
}
/* }}} */
/* plot method {{{ */
void GPE1D::plot(int nmodes, std::string &name) {
    std::ofstream file("/tmp/psi_x.txt");
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
        for(int i=0;i<_n;i++) {
            file << i*_dx-_xmax << ' ' << _psi[i].real() << ' ' << _psi[i].imag();
            for(int j=0;j<nmodes;j++) {
                file << ' ' << u[j*_n+i].real()
                    << ' ' << v[j*_n+i].real();
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi_x.txt\" !" << std::endl;
    }
    //Fourier space output
    std::ofstream file2("/tmp/psi_p.txt");
    if(file2.is_open()) {
        fftw_execute(_planFFT);
        const std::complex<double> *psi=_psip;
        for(int i=0;i<_n;i++) {
            double k=0.5/_xmax*(i<_n/2?i:i-_n);
            file2 << k << ' ' << psi[i].real() << ' ' << psi[i].imag();
            file2 << '\n';
        }
        file2.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psi_p.txt\" !" << std::endl;
    }
    //Gnuplot interface
    std::cout << "set size 1,1;set origin 0,0;unset key;"
        << "set multiplot layout 1,2;"
        << "set xlabel \"x\";set ylabel \"Density\";"
        << "plot \"/tmp/psi_x.txt\" using 1:($2**2+$3**2);"
        << "set xlabel \"k_x\";set ylabel \"Density\";"
        << "plot \"/tmp/psi_p.txt\" using 1:($2**2+$3**2);"
        << "unset multiplot;";
    /*
    for(int i=0;i<nmodes;i++) {
        std::cout << ",\"\" using 1:(2*$2*($" << 2*i+4 << "+$" << 2*i+5 << ")) title \"\"";
    }
    std::cout << ";";
    */
    std::cout << "pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* setHeader method {{{ */
void GPE1D::setHeader(std::ofstream &file) const {
    file.write((const char*)_type,3*sizeof(char));
    file.write((const char*)&_n,sizeof(int));
    file.write((const char*)&_xmax,sizeof(double));
    file.write((const char*)&_dx,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state GPE1D::getHeader(std::ifstream &file) {
    if(GPE::getHeader(file)==error)
        return error;
    int n;
    file.read((char*)&n,sizeof(int));
    file.read((char*)&_xmax,sizeof(double));
    file.read((char*)&_dx,sizeof(double));
    file.read((char*)&_mu,sizeof(double));
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
        double kx=2*pi/(_dx*_n)*(i<_n/2?i:i-_n);
        double k2=-kx*kx;
        _phase[i]=scale*exp(dt*_kterm*k2);
    }
}
/* }}} */
/* initialize method {{{ */
void GPE1D::initialize(Expression *pot) {
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi);
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    _planFFT=fftw_plan_dft_3d(1,1,_n,rspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    _planIFFT=fftw_plan_dft_3d(1,1,_n,pspace,rspace,FFTW_BACKWARD,FFTW_MEASURE);
    //Diagonal part
    VarDef vars;
    vars["X"]=new Constant(0);
    double *v=new double[_n];
    for(int i=0;i<_n;i++) {
        double x=i*_dx-_xmax;
        vars["X"]->set(&x);
        double vpot=*((double*)(pot->evaluate(vars)));
        vpot*=_vterm;
        _vpot[i]=vpot;
        v[i]=_kterm*(-2/(_dx*_dx))+vpot;
        _psi[i]=vpot<1?sqrt(1-vpot):0;
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
}
/* }}} */
/* ekin method {{{ */
double GPE1D::ekin() {
    fftw_execute(_planFFT);
    double res=0.;
    double n2=0.;
    for(int i=0;i<_n;i++) {
        double kx=2*pi/(_dx*_n)*(i<_n/2?i:i-_n);
        double psi2=std::norm(_psip[i]);
        res+=-_kterm*kx*kx*psi2;
        n2+=psi2;
    }
    return res/n2;
}
/* }}} */
/* }}} */
/* gpe.cpp */
