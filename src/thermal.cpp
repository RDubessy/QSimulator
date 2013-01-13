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
#include <iostream>
#include <common.h>     //For ConfigMap.
#include <expression.h> //For Expression, VarDef.
#include "thermal.h"
Thermal::Thermal(ConfigMap &config, VarDef &params) {
    _Ntot=getConfig(config,std::string("thermal::N"),1.0);
    _Nbec=getConfig(config,std::string("thermal::Nbec"),_Ntot);
    Expression *beta=parseString(config["thermal::beta"]);
    Expression *lambda=parseString(config["thermal::lambda"]);
    _beta=*((double*)(beta->evaluate(params)));
    _lambda=*((double*)(lambda->evaluate(params)));
    std::cerr << "[I] Adding thermal component:" 
        << "\n[I]\tbeta  =" << _beta
        << "\n[I]\tlambda=" << _lambda
        << std::endl;
    _n0=0;
    _Ntherm=0;
    //_Nbec=_Ntot;
}
Thermal::~Thermal() {
    if(_n0!=0)
        delete _n0;
}
/* thermal.cpp */
