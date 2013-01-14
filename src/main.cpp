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
/*!\todo Generalization to cylindrical coordinates.
 */
#include <iostream>     //For cerr/cout/endl...
#include <fstream>      //For ofstream/ifstream...
#include <common.h>     //For ConfigMap.
#include <expression.h> //For Expression, VarDef.
#include "thermal.h"    //For Thermal class.
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
    if(config.find("thermal::N")!=config.end()) {
        if(pot->find("X")) {
            if(pot->find("Y")) {
                gpe=new GPE2DThermal(config,eqn,pot,params);
            }
        } else if(pot->find("R"))
            gpe=new Polar1DThermal(config,eqn,pot,params);
    } else {
        if(pot->find("X")) {
            if(pot->find("Y")) {
                if(pot->find("Z")) {
                    if(eqn->find("LZ"))
                        gpe=new GPE3DROT(config,eqn,pot);
                    else
                        gpe=new GPE3D(config,eqn,pot);
                } else {
                    if(eqn->find("LZ"))
                        gpe=new GPE2DROT(config,eqn,pot);
                    else
                        gpe=new GPE2D(config,eqn,pot);
                }
            } else
                gpe=new GPE1D(config,eqn,pot);
        } else if(pot->find("R"))
            gpe=new Polar1D(config,eqn,pot);
    }
    if(gpe==0) {
        cerr << "[E] Unknown problem type!" << std::endl;
        return -1;
    }
    gpe->initialize(pot);
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
    if(config.find("imprint")!=config.end()) {
        int l=getConfig(config,string("imprint"),1);
        std::cerr << "[I] Imprinting a l=" << l << " circulation" << endl;
        gpe->imprint(l);
    }
    int verbose=0;
    if(config.find("verb")!=config.end()) {
        std::cerr << "[I] Verbose output" << std::endl;
        verbose=1;
    }
    if(config.find("groundstate")!=config.end()) {
        double dt=getConfig(config,string("general::dt"),1e-3);
        double tol=getConfig(config,string("general::tol"),1e-12);
        double dttol=getConfig(config,string("general::dttol"),0.999);
        gpe->findGroundState(dt,tol,dttol,log,verbose);
        gpe->save(out);
    }
    if(config.find("spectrum")!=config.end()) {
        range<int> def={0,0,1};
        range<int> m=getConfig(config,string("spectrum"),def);
        for(int i=m.min;i<=m.max;i+=m.incr)
            gpe->spectrum(log,i);
    }
    if(config.find("evolve")!=config.end()) {
        range<double> def={0.,1.,1e-3};
        range<double> t=getConfig(config,string("evolve"),def);
        gpe->evolve(t.min,t.incr,t.max,out);
    }
    if(config.find("plot")!=config.end()) {
        int n=5; //getConfig(config,string("plot"),0);
        gpe->plot(n,config["plot"]);
    }
    if(config.find("measure")!=config.end()) {
        std::cout << gpe->measure() << std::endl;
    }
    return 0;
};
/* }}} */
/* usage method {{{ */
void usage(const char *name) {
    cerr << "Usage: " << name << " CONFIG_FILE [OPTIONS] [ACTIONS]\n"
        << "The << qsimu >> program is intended to study properties of the Non-Linear\n"
        << "Schrodinger Equation (NLSE) also known in the cold atom community as the\n"
        << "Gross-Pitaevskii Equation (GPE). As of now it can compute the groundstate\n"
        << "of the equation, its time dynamics evolution and the spectrum of linear\n"
        << "(Bogoliubov) excitations.\n"
        << endl;
    cerr << "Possible [OPTIONS] are:\n"
        << "  --usage,--help    Display this screen and exits.\n"
        << "  --in=FILE         Load the initial state from file FILE.\n"
        << "  --out=FILE        Save the final state to file FILE.\n"
        << "  --log[=FILE]      Store the logs in the file log.txt or in the provided\n"
        << "                    (optionnal) FILE.\n"
        << "  --FIELD::KEYWORD=VALUE\n"
        << "                    Special option format which allows you to modify options\n"
        << "                    passed in CONFIG_FILE, namely assign a new VALUE to\n"
        << "                    KEYWORD in part FIELD of the file. See the documentation\n"
        << "                    for more details.\n"
        << endl;
    cerr << "Possible [ACTIONS] are:\n"
        << "  --groundstate     Computes the groundstate of the system and store the\n"
        << "                    result accordingly to --out=FILE. Use --log to get info\n"
        << "                    on the convergence towards the final state.\n"
        << "  --spectrum=RANGE\n"
        << "  --evolve=RANGE\n"
        << "  --plot\n"
        << "  --measure\n"
        << endl;
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
