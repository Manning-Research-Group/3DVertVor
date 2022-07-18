// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef POLYGON_H_INCLUDED
#define POLYGON_H_INCLUDED

class Polygon;
#include "../Run/Run.h"
#include "../Edge/Edge.h"
#include "../Cell/Cell.h"

class Polygon {
public:
    long int id_;
    double center_[3];
    double area_;
    double tension_;
    double volumeForce_[3];
    double interfaceForce_[3];
    double cellcenter_[3];
    double normvector_[3];
    double sigma_;
    int twincell_[4];
    int type_;
    std::vector<Edge *> edges_;
    std::vector<Vertex *> vertices_;
    std::vector<Cell *> cells_;
    explicit Polygon(Run *, long int);

    int updateVertices();
    int updateCenter();
    int updateArea();
    bool crossBoundary();
    bool checkH();
    int shrink(Edge *);
    int expand(Edge *);
    int logEdges(std::string);
    int SpringGon(double curtime);
private:
    Run * run_;
};

#endif
