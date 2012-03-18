/*!\mainpage Quantum simulation software.
 *
 * \todo Generalization to higher dimensions (cartesian or cylindrical/spherical).
 */
#include <iostream>     //For cerr/cout/endl...
#include <fstream>      //For ofstream/ifstream...
#include <common.h>     //For ConfigMap.
#include <expression.h> //For Expression, VarDef.
#include <cvm.h>        //For cvm matrices.
#include "gpe.h"        //For GPE class (and sub-classes).
/* mainFunction method {{{ */
int mainFunction(ConfigMap &config) {
    if(config.find("general::equation")==config.end()) {
        cerr << "[E] No problem defined!" << endl;
        return -1;
    }
    cerr << "[I] Initializing problem...\n";
    Expression *eqn=parseString(config["general::equation"]);
    Expression *pot=parseString(config["general::potential"]);
    cerr << "[I] Equation: " << eqn << endl;
    cerr << "[I] Potential: " << pot << endl;
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
    GPE *gpe=0;
    if(config["general::type"]=="Polar") {
        gpe=new Polar1D(config,eqn,pot);
    } else if(config["general::type"]=="Cartesian1D") {
        gpe=new GPE1D(config,eqn,pot);
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
    string log="";
    if(config.find("log")!=config.end()) {
        log="log.txt";
        if(config["log"].size()>0) {
            log=config["log"];
        }
    }
    if(config.find("groundstate")!=config.end()) {
        double dt=getConfig(config,string("general::dt"),1e-3);
        double tol=getConfig(config,string("general::tol"),1e-8);
        double dttol=getConfig(config,string("general::dttol"),0.9999);
        gpe->findGroundState(dt,tol,dttol,log);
        gpe->save(out);
    }
    if(config.find("plot")!=config.end()) {
        int n=getConfig(config,string("plot"),0);
        gpe->plot(n);
    }
    if(config.find("spectrum")!=config.end()) {
        range<int> def={0,0,1};
        range<int> m=getConfig(config,string("spectrum"),def);
        for(int i=m.min;i<=m.max;i+=m.incr)
            gpe->spectrum(log,i);
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
    if(config.find("usage")!=config.end()||config.find("help")!=config.end()) {
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
