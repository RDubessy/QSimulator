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
#include "gpe.h"
#include <fstream>      //For ifstream, ofstream
#include <iomanip>      //For setfill, setw
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
    _mu=0;
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
 * \param verb an integer controlling the verbosity of the displayed output.
 */
void GPE::findGroundState(double dttest, double tol, double dttol, string &name, int verb) {
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
                double deltamu=-std::log(tmp)/dt;
                _mu+=deltamu;
                eps=std::abs(deltamu/_mu);
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
    std::cout.precision(15);
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
    vars["LZ"]=(Constant*)&zero;
    vars["t_"]=new Constant(t);
    //Kinetic term pre-factor
    vars["DELTA"]=(Constant*)&one;
    vars["VEXT"]=(Constant*)&zero;
    vars["RHO"]=(Constant*)&zero;
    _kterm=*((double*)(_H->evaluate(vars)));
    //Interaction term pre-factor
    vars["DELTA"]=(Constant*)&zero;
    vars["VEXT"]=(Constant*)&zero;
    vars["RHO"]=(Constant*)&one;
    _gN=*((double*)(_H->evaluate(vars)));
    //Potential term pre-factor
    vars["DELTA"]=(Constant*)&zero;
    vars["VEXT"]=(Constant*)&one;
    vars["RHO"]=(Constant*)&zero;
    _vterm=*((double*)(_H->evaluate(vars)));
}

/* }}} */
/* measure method {{{ */
std::string GPE::measure() {
    std::ostringstream mes;
    mes.precision(15);
    mes << _mu << ' ' << epot() << ' ' << ekin();
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
/* getHeader method {{{ */
state GPE::getHeader(std::ifstream &file) {
    char type[3];
    file.read((char*)type,3*sizeof(char));
    if(type[0]!=_type[0] || type[1]!=_type[1] || type[2]!=_type[2]) {
        std::cerr << "[E] Incompatible type !" << std::endl;
        return error;
    }
    return ok;
}
/* }}} */
/* }}} */
/* gpe.cpp */
