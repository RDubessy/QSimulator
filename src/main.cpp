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
    if(config["general::equation"].size()==0) {
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
    if(config["groundstate"].size()>0) {
        double dt=getConfig(config,string("general::dt"),1e-3);
        double tol=getConfig(config,string("general::tol"),1e-8);
        gpe->findGroundState(dt,tol);
        gpe->save(out);
    }
    if(config["plot"].size()>0) {
        int n=getConfig(config,string("plot"),0);
        gpe->plot(n);
    }
    if(config["spectrum"].size()>0) {
        int m=getConfig(config,string("spectrum"),0);
        for(int i=0;i<=m;i++)
            gpe->spectrum(i);
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
