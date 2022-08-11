// Copyright June 2021 Tao Zhang @ Shanghai Jiao Tong University.  All Rights Reserved.
// Author: Tao Zhang @ Shanghai Jiao Tong University, zhangtao.scholar@sjtu.edu.cn
// Corresponding author: Jennifer Schwarz @ Syracuse University, jschwarz@physics.syr.edu

#ifndef RUN_H_INCLUDED
#define RUN_H_INCLUDED

class Run;
#include <vector>
#include "../Vertex/Vertex.h"
#include "../Edge/Edge.h"
#include "../Polygon/Polygon.h"
#include "../Cell/Cell.h"
#include "../Energy/Volume.h"
#include "../Energy/Interface.h"
#include "../Reconnection/Reconnection.h"

#include "misc/geometry/Vector3D.h"
#include "misc/geometry/Matrix3x3.h"
#include "misc/geometry/Matrix3x3x3.h"
#include "misc/geometry/Ellipsoid.h"
#include "misc/geometry/EllipsoidByUnitPointMassPolyhedron.h"

class Run {
  public:
    double  dt_;    // integration time step
    double  dtr_;   // time interval of network reconnection
    double  mu_;   // inverse damping coefficient of vertex
    double  kB_;
    double  temperature_;
    double  temperaturebot_;
    double  Lx_;
    double  Ly_;
    double  Lz_;
    int     NCell_;
    double  simulation_time_;
    double   t_start_;
    double   t_end_;
    long int count_dump_;
    double   dump_period_;
    long int count_log_;
    double   log_period_;
    long int count_reconnect_;
    long int count_vertices_;
    long int count_edges_;
    long int count_polygons_;
    long int count_cells_;
    double sigma_;


    double s01_;
    double s02_ = 6.70;
    double s03_ = 6.70;
    double s04_;
    double s05_;
    double s06_;
    double s07_;
    double v01_;
    double v02_ = 1.3;
    double v03_ = 1.3;
    double v04_;
    double v05_;
    double v06_;
    double v07_;
    double l0s_;

    double sigma12_;
    double sigma23_ = 0.5;
    double sigma31_ = 3.;
    double sigma14_;
    double sigma24_ = 0.0;
    double sigma34_ = 0.0;
    double sigma15_;
    double sigma25_ = 0.0;
    double sigma35_ = 0.0;
    double sigma45_;
    double sigma16_ = 0.0;
    double sigma26_ = 0.5;
    double sigma36_ = 2.;
    double sigma46_;
    double sigma56_;
    double sigma17_;
    double sigma27_ = 0.0;
    double sigma37_ = 0.0;
    double sigma47_;
    double sigma57_;
    double sigma67_;

    int placodeon_ = 0;
    double placodex_ = 0.0;
    double placodey_ = 0.0;
    int placodesize_ = 0;

    int type1count;
    int type2count;

    Volume * volume_;
    Interface * interface_;
    Reconnection * reconnection_;
    std::stringstream verboseReconnection_;

    std::vector<Vertex *> vertices_;
    std::vector<Edge *> edges_;
    std::vector<Polygon *> polygons_;
    std::vector<Cell *> cells_;

    Run();
    int     start();
    int     updatePolygonVertices();
    int     updatePolygonCells();
    int     updateCellVertices();
    int     updateCellShapeIndex();
    int     updateVertexEdges();
    int     updateVertexCells();
    int     updateGeoinfo();
    int     updateVerticesVelocity();
    int     updateVerticesPosition();
    int     deleteVertex(Vertex *);
    int     deleteEdge(Edge *);
    int     deletePolygon(Polygon *);
    int     resetPosition(double *);
    Edge *  addEdge(Vertex *, Vertex *);
    int     dumpConfigurationVtk();
    int     dumpCellCenter();
    int     dumpCellShapeIndex();
    int     dumpTopo();
    int     dumpReconnection();
    int     dumpDemix();
    int     InitializeCellType();
    int     InitializeCellUniform();
    int     InitializeCellDelam();
    int     InitializePlacode();
    int     PlacodeCenter();
    int     InitalizeSolidifcation();
    int     InitializeCellDirectors();
    int     InitalizeErin();
    double  Randnormal(double stdev);
    EllipsoidByUnitPointMassPolyhedron fitEllipsoidByUnitPointMassPolyhedron(int cellID);
    Matrix3x3 _unitPointMassMomentOfInertiaTensor;
};

#endif