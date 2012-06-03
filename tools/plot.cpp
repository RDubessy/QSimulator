#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <complex>
using namespace std;
struct param {
    double mu;
    double dr;
    double rmin;
    double *psi;
    complex<double> *u;
    complex<double> *v;
    int n;
    int l;
    int m;
};
double pi=acos(-1);
void printUsage(const char *name) {
    cerr << "Usage:\n" << name << " file l optfile m\n"
        << "Where:\n"
        << "\tfile contains the wavefunction data\n"
        << "\tl is the circulation in the groundstate\n"
        << "\toptfile contains the excitation spectrum\n"
        << "\tm is the circulation of the excitation\n";
    cerr.flush();
};
void readPsi(const char *name, param &p) {
    ifstream file(name);
    if(!file.is_open() || p.psi!=0)
        return;
    /* get headers */
    char type[3];
    file.read((char*)type,3*sizeof(char));
    file.read((char*)&p.n,sizeof(int));
    file.read((char*)&p.rmin,sizeof(double));
    file.read((char*)&p.dr,sizeof(double));
    file.read((char*)&p.mu,sizeof(double));
    if(type[0]!='P' || type[1]!='o' || type[2]!='l') {
        cerr << "[E] Type: " << type[0] << type[1] << type[2] << " not supported!" << endl;
        return;
    }
    /* read data */
    cerr << "[I] Loading wavefunction (" << p.n << ")" << endl;
    p.psi=new double[p.n];
    for(int i=0;i<p.n;i++) {
        file.read((char*)&(p.psi[i]),sizeof(double));
    }
    file.close();
    return;
};
void readExc(const char *name, param &p) {
    ifstream file(name);
    if(!file.is_open() || p.psi!=0)
        return;
    /* get headers */
    char type[3];
    double rmin,dr,mu;
    file.read((char*)type,3*sizeof(char));
    file.read((char*)&p.n,sizeof(int));
    file.read((char*)&rmin,sizeof(double));
    file.read((char*)&dr,sizeof(double));
    file.read((char*)&mu,sizeof(double));
    if(type[0]!='P' || type[1]!='o' || type[2]!='l') {
        cerr << "[E] Type: " << type << " not supported !" << endl;
        return;
    }
    /* read data */
    p.u=new complex<double>[p.n*p.n];
    p.v=new complex<double>[p.n*p.n];
    for(int i=0;i<p.n;i++) {
        complex<double> e;
        file.read((char *)&e,sizeof(complex<double>));
        for(int j=0;j<p.n;j++)
            file.read((char *)&p.u[i*p.n+j],sizeof(complex<double>));
        for(int j=0;j<p.n;j++)
            file.read((char *)&p.v[i*p.n+j],sizeof(complex<double>));
    }
    file.close();
    return;
};
void analyse(param &p) {
    //Average r value
    double r0=0.;
    for(int i=0;i<p.n;i++) {
        double r=p.rmin+i*p.dr;
        r*=p.psi[i];
        r0+=r*r;
    }
    r0*=2*pi*p.dr;
    cerr << "Groundstate properties:\n"
        << "\tChemical potential: " << p.mu << "\n"
        << "\tAverage radius <r>: " << r0 << "\n";
    cerr.flush();
    if(p.u!=0 && p.v!=0) {
        double r2=0;
        double r1=0;
        for(int i=0;i<p.n;i++) {
            double r=p.rmin+i*p.dr;
            double arg=p.psi[i]*(p.u[i].real()+p.v[i].real());
            r2+=r*r*arg;
            r1+=r*arg;
        }
        cerr << "Excited state properties:\n"
            << "\tAverage radius <r>: " << r2/r1 << "\n";
        cerr.flush();
    }
    return;
};
int main(int argc, char *argv[]) {
    param p={0.,0.,0.,0,0,0,0,0,0};
    switch(argc) {
        case 5:
            p.m=atoi(argv[4]);
        case 4:
            readExc(argv[3],p);
        case 3:
            p.l=atoi(argv[2]);
        case 2:
            readPsi(argv[1],p);
            break;
        case 1:
        default:
            printUsage(argv[0]);
            return -1;
    }
    if(p.psi==0)
        return 0;
    analyse(p);
    ofstream file("/tmp/psi2D.txt");
    int nj=abs(p.m)*20;
    if(nj<256) nj=256;
    if(p.u!=0 && p.v!=0) {
        for(int j=0;j<nj;j++) {
            for(int i=0;i<p.n;i++)
                file << p.rmin+i*p.dr << ' ' << j*2*pi/(nj-1) 
                    << ' ' << p.psi[i] 
                    << ' ' << p.u[i].real()
                    << ' ' << p.u[i].imag()
                    << ' ' << p.v[i].real()
                    << ' ' << p.v[i].imag()
                    << '\n';
            file << '\n';
        }
    } else {
        for(int j=0;j<nj;j++) {
            for(int i=0;i<p.n;i++)
                file << p.rmin+i*p.dr << ' ' << j*2*pi/(nj-1) 
                    << ' ' << p.psi[i] << '\n';
            file << '\n';
        }
    }
    file.close();
    cout << "set mapping cylindrical;set pm3d;unset surface;set size square;set view map;";
    if(p.u!=0&&p.v!=0) {
        cout << "unset key;"
            << "set multiplot layout 1,2;"
            << "splot \"/tmp/psi2D.txt\" using 2:(2*$3*($4*cos($2*(" << p.m << "))-$5*sin($2*(" << p.m << "))+$6*cos($2*(" << p.m << "))-$7*sin($2*(" << p.m << ")))):1;"
            << "splot \"/tmp/psi2D.txt\" using 2:(arg($3+($4+{0,1}*$5)*exp({0,1}*" << p.m << "*$2)+($6-{0,1}*$7)*exp({0,-1}*" << p.m << "*$2))+pi):1;"
            << "unset multiplot;";
    } else {
        cout << "splot \"/tmp/psi2D.txt\" using 2:($3*$3):1 title \"\";";
    }
    cout << "pause mouse" << endl;
    return 0;
};
/* plot.cpp */
