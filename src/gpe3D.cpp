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
#include "gpe3D.h"
#include <fstream>      //For ifstream, ofstream
/* class GPE3D implementation {{{ */
/* Constructor {{{ */
GPE3D::GPE3D(ConfigMap &config, Expression *H, Expression *pot) : GPE(H,pot) {
    std::cerr << "[I] Initializing a 3D Cartesian system" << std::endl;
    _type="x3D";
    _nx=getConfig(config,string("grid::nx"),64);
    _ny=getConfig(config,string("grid::ny"),_nx);
    _nz=getConfig(config,string("grid::nz"),_ny);
    _xmax=getConfig(config,string("grid::Lx"),1.);
    _ymax=getConfig(config,string("grid::Ly"),_xmax);
    _zmax=getConfig(config,string("grid::Lz"),_zmax);
    _dx=_xmax/(_nx-1);
    _dy=_ymax/(_ny-1);
    _dz=_zmax/(_nz-1);
    _xmax/=2;
    _ymax/=2;
    _zmax/=2;
    int n=_nx*_ny*_nz;
    _psi.resize(n);
    allocate(n);
}
/* }}} */
/* spectrum method {{{ */
void GPE3D::spectrum(string &name, int m) {
    std::cerr << "[E] spectrum method NOT implemented for a 3D system !" << std::endl;
    return;
}
/* }}} */
/* norm methods {{{ */
double GPE3D::norm(cvm::rvector &psi) const {
    return sqrt(_dx*_dy*_dz)*psi.norm2();
}
double GPE3D::norm(cvm::cvector &psi) const {
    return sqrt(_dx*_dy*_dz)*psi.norm2();
}
/* }}} */
/* plot method {{{ */
void GPE3D::plot(int nmode, std::string &name) {
    std::ofstream file;
    file.open("/tmp/psiXY.txt");
    if(file.is_open()) {
        const std::complex<double> *psi=_psi.get();
        for(int j=0;j<_ny;j++) {
            for(int i=0;i<_nx;i++) {
                double psi2=0;
                for(int k=0;k<_nz;k++)
                    psi2+=std::norm(psi[i+(j+k*_ny)*_nx]);
                psi2/=2*_zmax;
                file << i*_dx-_xmax << ' ' << j*_dy-_ymax << ' ' << psi2
                    << '\n';
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psiXY.txt\" !" << std::endl;
    }
    file.open("/tmp/psiXZ.txt");
    if(file.is_open()) {
        const std::complex<double> *psi=_psi.get();
        for(int k=0;k<_nz;k++) {
            for(int i=0;i<_nx;i++) {
                double psi2=0;
                for(int j=0;j<_ny;j++)
                    psi2+=std::norm(psi[i+(j+k*_ny)*_nx]);
                psi2/=2*_ymax;
                file << k*_dz-_zmax << ' ' << i*_dx-_xmax << ' ' << psi2
                    << '\n';
            }
            file << '\n';
        }
        file.close();
    } else {
        std::cerr << "[E] Can't open file \"/tmp/psiXZ.txt\" !" << std::endl;
    }
    std::cout << "set view map;unset surface;set pm3d;unset key;"
        << "set size square;set size 1,1;set origin 0,0;"
        << "set multiplot layout 1,2;"
        << "splot \"/tmp/psiXY.txt\" using 1:2:3;"
        << "splot \"/tmp/psiXZ.txt\" using 1:2:3;"
        << "unset multiplot;"
        << "pause mouse\n";
    std::cout.flush();
    return;
}
/* }}} */
/* setHeader method {{{ */
void GPE3D::setHeader(std::ofstream &file) const {
    file.write((const char*)_type,3*sizeof(char));
    file.write((const char*)&_nx,sizeof(int));
    file.write((const char*)&_ny,sizeof(int));
    file.write((const char*)&_nz,sizeof(int));
    file.write((const char*)&_xmax,sizeof(double));
    file.write((const char*)&_ymax,sizeof(double));
    file.write((const char*)&_zmax,sizeof(double));
    file.write((const char*)&_dx,sizeof(double));
    file.write((const char*)&_dy,sizeof(double));
    file.write((const char*)&_dz,sizeof(double));
    file.write((const char*)&_mu,sizeof(double));
}
/* }}} */
/* getHeader method {{{ */
state GPE3D::getHeader(std::ifstream &file) {
    if(GPE::getHeader(file)==error)
        return error;
    int nx,ny,nz;
    file.read((char*)&nx,sizeof(int));
    file.read((char*)&ny,sizeof(int));
    file.read((char*)&nz,sizeof(int));
    file.read((char*)&_xmax,sizeof(double));
    file.read((char*)&_ymax,sizeof(double));
    file.read((char*)&_zmax,sizeof(double));
    file.read((char*)&_dx,sizeof(double));
    file.read((char*)&_dy,sizeof(double));
    file.read((char*)&_dz,sizeof(double));
    file.read((char*)&_mu,sizeof(double));
    if(nx!=_nx || ny!=_ny || nz!=_nz) {
        std::cerr << "[E] Incompatible size !" << std::endl;
        return error;
    }
    return ok;
}
/* }}} */
/* computePhase method {{{ */
void GPE3D::computePhase(std::complex<double> dt) {
    double scale=1./(_nx*_ny*_nz);
    for(int k=0;k<_nz;k++) {
        double kz=2*pi/(_dz*_nz)*(k<_nz/2?k:k-_nz);
        double k2z=kz*kz;
        for(int j=0;j<_ny;j++) {
            double ky=2*pi/(_dy*_ny)*(j<_ny/2?j:j-_ny);
            double k2y=ky*ky;
            for(int i=0;i<_nx;i++) {
                double kx=2*pi/(_dx*_nx)*(i<_nx/2?i:i-_nx);
                double k2x=kx*kx;
                double k2=-(k2x+k2y+k2z);
                _phase[i+(j+k*_ny)*_nx]=scale*std::exp(dt*_kterm*k2);
            }
        }
    }
}
/* }}} */
/* initialize method {{{ */
void GPE3D::initialize(Expression *pot) {
    //Only a diagonal part
    VarDef vars;
    vars["X"]=new Constant(0);
    vars["Y"]=new Constant(0);
    vars["Z"]=new Constant(0);
    int n=_nx*_ny*_nz;
    double *psi=new double[n];
    for(int k=0;k<_nz;k++) {
        double z=k*_dz-_zmax;
        vars["Z"]->set(&z);
        for(int j=0;j<_ny;j++) {
            double y=j*_dy-_ymax;
            vars["Y"]->set(&y);
            for(int i=0;i<_nx;i++) {
                double x=i*_dx-_xmax;
                vars["X"]->set(&x);
                double vpot=*((double*)(pot->evaluate(vars)));
                vpot*=_vterm;
                _vpot[i+(j+k*_ny)*_nx]=vpot;
                psi[i+(j+k*_ny)*_nx]=vpot<1?sqrt(1-vpot):0;
            }
        }
    }
    initializeFFT();
    _psi=cvm::cvector(psi,n);
    delete[] psi;
}
/* }}} */
/* initializeFFT method {{{ */
void GPE3D::initializeFFT() {
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi.get());
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    _planFFT=fftw_plan_dft_3d(_nz,_ny,_nx,rspace,pspace,FFTW_FORWARD,FFTW_MEASURE);
    _planIFFT=fftw_plan_dft_3d(_nz,_ny,_nx,pspace,rspace,FFTW_BACKWARD,FFTW_MEASURE);
}
/* }}} */
/* ekin method {{{ */
double GPE3D::ekin() {
    double res=0.;
    double n2=0.;
    for(int k=0;k<_nz;k++) {
        double k2z=(cos((2*pi*k)/_nz)-1)/(_dz*_dz);
        for(int j=0;j<_ny;j++) {
            double k2y=(cos((2*pi*j)/_ny)-1)/(_dy*_dy);
            for(int i=0;i<_nx;i++) {
                double k2x=(cos((2*pi*i)/_nx)-1)/(_dx*_dx);
                double psi2=std::norm(_psip[i+(j+k*_ny)*_nx]);
                res+=2*_kterm*(k2x+k2y+k2z)*psi2;
                n2+=psi2;
            }
        }
    }
    return res/n2;
}
/* }}} */
/* }}} */
/* class GPE3DROT implementation {{{ */
/* Constructor {{{ */
GPE3DROT::GPE3DROT(ConfigMap &config, Expression *H, Expression *pot) : GPE3D(config,H,pot) {
    _phase2=new std::complex<double>[_nx*_ny*_nz];
    update(0);
}
/* }}} */
/* initializeFFT method {{{ */
void GPE3DROT::initializeFFT() {
    std::cerr << "[I] FFT in Rotating Frame" << std::endl;
    //Initialize the plans for fast fourier transform
    fftw_complex *rspace=reinterpret_cast<fftw_complex*>(_psi.get());
    fftw_complex *pspace=reinterpret_cast<fftw_complex*>(_psip);
    fftw_iodim dimx,dimy,dimz;
    dimz.n=_nz;
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
void GPE3DROT::computePhase(std::complex<double> dt) {
    double scale=1./(_nx*_ny*_nz);
    double xo=0.5*(_nx-1);
    double yo=0.5*(_ny-1);
    for(int k=0;k<_nz;k++) {
        double kz=2*pi/(_dz*_nz)*(k<_nz/2?k:k-_nz);
        double k2z=-kz*kz;
        for(int j=0;j<_ny;j++) {
            double y=j-yo;
            double ky=2*pi/(_dy*_ny)*(j<_ny/2?j:j-_ny);
            double k2y=-ky*ky;
            for(int i=0;i<_nx;i++) {
                double x=i-xo;
                double kx=2*pi/(_dx*_nx)*(i<_nx/2?i:i-_nx);
                double k2x=-kx*kx;
                _phase[i+(j+k*_ny)*_nx]=scale*std::exp(dt*(_kterm*(k2x+k2z)-_oterm*y*kx));
                _phase2[i+(j+k*_ny)*_nx]=std::exp(dt*(_kterm*k2y+_oterm*x*ky));
            }
        }
    }
    return;
}
/* }}} */
/* doStep method {{{ */
void GPE3DROT::doStep(std::complex<double> dt) {
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
void GPE3DROT::update(double t) {
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
std::string GPE3DROT::measure() {
    //Compute <Lz> & rotation energy
    double lz1=0.;
    double lz2=0.;
    double xo=0.5*(_nx-1);
    double yo=0.5*(_ny-1);
    double n2=0.;
    //Term ykx
    fftw_execute(_planFFTxz);
    for(int k=0;k<_nz;k++) {
        for(int j=0;j<_ny;j++) {
            double y=j-yo;
            for(int i=0;i<_nx;i++) {
                double kx=sin((2*pi*i)/_nx);
                double psi2=std::norm(_psip[i+(j+k*_ny)*_nx]);
                lz1-=y*kx*psi2;
                n2+=psi2;
            }
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
    for(int k=0;k<_nz;k++) {
        for(int j=0;j<_ny;j++) {
            double ky=sin((2*pi*j)/_ny);
            for(int i=0;i<_nx;i++) {
                double x=i-xo;
                double psi2=std::norm(_psip[i+(j+k*_ny)*_nx]);
                lz2+=x*ky*psi2;
                n2+=psi2;
            }
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
