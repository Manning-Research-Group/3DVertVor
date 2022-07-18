// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef INTERFACE_H_INCLUDED
#define INTERFACE_H_INCLUDED

class Interface;
#include "../Run/Run.h"

class Interface {
public:
    double s0_;

    double energy_;

    explicit Interface(Run *);

    int updateForces();
    int updatePolygonForces(Polygon *);
    int updateTension();
    int updateEnergy();
private:
    Run * run_;
};

#endif
