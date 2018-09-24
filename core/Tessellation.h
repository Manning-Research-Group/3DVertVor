#ifndef TESSELLATION_H
#define	TESSELLATION_H

#include <vector>
#include <string>
#include <unordered_map>
#include <functional>

#if USE_NETCDF
#include <netcdfcpp.h>
#endif

#include "misc/other/Random.h"

#include "misc/geometry/Vector3D.h"
#include "misc/geometry/SphericalHarmonicsForVector.h"

#include "misc/minimization/MinimizerWithDerivative.h"
#include "dynamicalMatrix/DynamicalMatrixOld.h"
#include "dynamicalMatrix/DynamicalMatrixTest.h"
#include "povray/PovrayGenerator.h"

#include "CellType.h"
#include "Cell.h"
#include "CellPositionMailbox.h"


class DirectedFace;
class VertexOfCell;
struct PeriodicBox;

class Tessellation {
public:
  Tessellation(PeriodicBox &box);
  virtual ~Tessellation();
  void cleanup();

  // queries
  const PeriodicBox &box() const { return _box; }
  const std::vector<Cell*> &cells() const { return _cells; }
  
  // creation
  void addCell(const CellType &cellType, const Vector3D &position);
  void addCellsAtRandomPositions(const CellType &cellType, const int NumberOfCells);
  
  // minimization
  void setDofs(MinimizerWithDerivative &minimizer);
  bool relax(MinimizerWithDerivative &minimizer);

  // dynamics
  double time() const { return _time; }
  void timeStep(const double deltaT);

  // topology and geometry
  bool topologyChanged();
  bool topologyChangedBruteForce();
  bool checkVoronoiTesselationFromVoroLibrary() const;
  virtual void computeGeometry();
  /** returns true if topology changed. */
  bool computeTopologyAndGeometry();
  void affinelyRescaleBox(double factorX, double factorY, double factorZ);
  Vector3D averagePositionWithoutBox() const {
    Vector3D sumPosition(0,0,0);
    for(const Cell *c : cells()) sumPosition += c->positionWithoutBox();
    return sumPosition/cells().size();
  }
  template<int L> double bondOrientationalOrder() const;
  void setEdgeRestLengthsTo(double l);
  void setEdgeRestLengthsToCurrentLengths();

  
  // energy
  virtual double energy() const;
  
  // forces
  virtual void computeEnergyDerivatives();
  virtual double stressXx() const;
  virtual double stressYy() const;
  virtual double stressZz() const;
  virtual double stressYx() const; // first: normal vector, second: force vector
  virtual double pressure() const;
  virtual double pressureAlt() const;
  Matrix3x3 numericalDerVertexPositionWrtCellPosition(VertexOfCell *v, Cell *c);
  Matrix3x3 numericalDerFaceAreaWrtCellPosition(DirectedFace *f, Cell *c);
  Vector3D numericalDerCellEnergyWrtCellPosition(Cell *ce, Cell *cp);
  Vector3D numericalForceOnCell(Cell *c);
  double totalCellForceNormSq() const;
  
  // dynamical matrix
  void computeAndSolveDynamicalMatrix();
  DynamicalMatrixOld *dynamicalMatrix() const { return _dynamicalMatrix; }
//  DynamicalMatrixTest *dynamicalMatrix() const { return _dynamicalMatrix; }
  double pureShearModulusXy(double eigenvalueCutoff) const { return 0.25*_dynamicalMatrix->modulusWithRespectToControlParameter(0, eigenvalueCutoff)/_box.volume(); }
  double simpleShearModulusYx(double eigenvalueCutoff) const { return _dynamicalMatrix->modulusWithRespectToControlParameter(1, eigenvalueCutoff)/_box.volume(); }
  double bulkModulus(double eigenvalueCutoff) const { return _dynamicalMatrix->modulusWithRespectToControlParameter(2, eigenvalueCutoff)*_box.volume(); }
  double cellPairModulus(const Cell *c1, const Cell *c2, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const;
//  double cellPairModulusNumerical(const Cell *c1, const Cell *c2);
  double cellPairModulusX(const Cell *c1, const Cell *c2, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const;
  double cellPairModulusY(const Cell *c1, const Cell *c2, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const;
  double cellPairModulusZ(const Cell *c1, const Cell *c2, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const;
  
  // drawing
  virtual void saveAsImageDefault(const std::string &filename, int width=800, int height=800, 
                                    const Vector3D &CameraPositionOffset=Vector3D(1.6, -0.8, 1.37), const Vector3D &CameraLooksAt=Vector3D(0.5, 0.5, 0.375),
                                    const Color &Background=Color(0.85,0.90,1.00)) const;
  void saveAsImage(const std::string &filename, std::function<void(PovrayGenerator &g)> drawElements, int width=800, int height=800,
                                    const Vector3D &CameraPositionOffset=Vector3D(1.6, -0.8, 1.37), const Vector3D &CameraLooksAt=Vector3D(0.5, 0.5, 0.375),
                                    const Color &Background=Color(0.85,0.90,1.00)) const;
  void drawCellPositionsDefault(PovrayGenerator &g, const Color &color=Color(0.7, 0.7, 0.7), double radius=0.4) const;
  void drawCells(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c)> drawCell) const;
  void drawEdgesDefault(PovrayGenerator &g, const Color &color=Color(1,0.4,0.45), double radius=0.05) const;
  void drawEdges(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c, const DirectedEdgeOfCell *e)> drawEdge) const;
  void drawVerticesDefault(PovrayGenerator &g, const Color &color=Color(0.8,1.0,0.3), double radius=0.06) const;
  void drawVertices(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c, const VertexOfCell *v)> drawVertex) const;

  // loading and saving; compatible with the LiuJamming code
  bool loadFromDump(const CellType &cellParameters, const std::string &filename);
  bool writeToDump(const std::string &filename);
  
#if USE_NETCDF
  /** this uses a format similar to the CStaticDatabase class of the LiuJamming code */
  bool loadFromNetCdfFileWithCellParameters(const std::string &Path, CellType &parameters, int record=0);
  /** this uses a format similar to the CStaticDatabase class of the LiuJamming code */
  bool saveAsNetCdfFile(const std::string &Path) const;
#endif
  
protected:
  std::vector<Cell*> _cells;

  // geometry
  CellPositionMailbox _mailbox;
  PeriodicBox &box() { return _box; }

  // dynamical matrix
  DynamicalMatrixOld *_dynamicalMatrix;
  Vector3D _boxDimensionsVector;
//  DynamicalMatrixTest *_dynamicalMatrix;
  virtual void computeDynamicalMatrix();
  void addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(DynamicalMatrixOld &dm, const Cell *c1, const Cell *c2, const Matrix3x3 &derivative) const;
  void addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(DynamicalMatrixOld &dm, int cp, const Cell *c, const Vector3D &derivative) const;

  // loading and saving; compatible with the LiuJamming code
  bool loadFromDumpHelper(const std::string &filename, std::function<Cell*(const Vector3D &initialPosition)> createCell);

private:
  const static double DistanceSqCutoffForPeriodicityInference;
  const static double NumericalDerivativeDifference;

  PeriodicBox &_box;
  double _time;

  // topology and geometry
  bool _topologicalElementsPresent;
  void resetTopology();
  void setTopologyByVoronoiTesselation();
  void computeVertexPositions();
  
  // dynamical matrix
  std::unordered_map<const Cell*,int> _cellToDynamicalMatrixIndex;
  void addCellEnergyDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const Matrix3x3 &derivative);
  friend class Cell;
  friend class DirectedEdgeOfCell;
};

#endif	/* TESSELATION_H */

