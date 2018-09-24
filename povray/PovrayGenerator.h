/* 
 * File:   PovrayGenerator.h
 * Author: arbeit
 *
 * Created on October 9, 2015, 9:35 PM
 */

#ifndef POVRAYGENERATOR_H
#define	POVRAYGENERATOR_H

#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include "misc/geometry/Vector3D.h"
#include "Color.h"
#include "PovrayMesh.h"

class PovrayGenerator {
public:
  PovrayGenerator();
  ~PovrayGenerator();
  void start(const Color &background=Color(0.85,0.90,1.00));
//  void start(const Color &background=Color(1,1,1));
  void create(std::string filename, int width=1024);

  void addLight(const Vector3D &position, const Color &color);
  void setCamera(const Vector3D &position, const Vector3D &looksAt, const double aspectRatio=4.0/3.0);
  void setVariable(const std::string &name, double value);

  // TODO: this is not pretty, maybe one should use RAII or so for this kind of stuff
  bool startUnion();
  bool startUnionWithBoxIntersection(const Vector3D &Corner1, const Vector3D &Corner2);
  void endUnion();
  void endUnionWithColor(const Color &color);
  void endUnionWithColorF(const Color &color, const double filterTransparency);
  void endUnionWithColorT(const Color &color, const double transmTransparency);

  void drawCoordinateSystem();
  void drawSphere(const Vector3D &c, const double radius);
  void drawCylinder(const Vector3D &p1, const Vector3D &p2, const double radius);
  void drawArrow(const Vector3D &p1, const Vector3D &p2, const double radius);

  template <typename T> void drawMesh(const PovrayMesh<T> &mesh) {
    _output << "  mesh2{\n";

    auto positions = mesh.positions();
    _output << "    vertex_vectors{ " << positions.size();
    for(const Vector3D &p : positions) {
      _output << ", <" << p << ">";
    }
    _output << "}\n";

    auto triangles = mesh.triangles();
    _output << "    face_indices{ " << triangles.size();
    for(const IntegerTriple &t : triangles) {
      _output << ", <" << t << ">";
    }
    _output << "}\n";

    _output << "    inside_vector <0,0,1>\n";
    _output << "  }\n";
  }
  template <typename T> bool startUnionWithMeshIntersection(const PovrayMesh<T> &mesh) {
    if(_currentlyInUnion) {
      return false;
    }

    _output << "union{intersection{\n";
    drawMesh(mesh);
    _output << "  union{\n";
    _currentlyInUnion = true;
    _currentlyInIntersectionUnion = true;
    return true;
  }
  
private:
  static const std::string TempFolderName;
  static const std::string TempFileNameMask;
  static const std::string LogFileName;
  static const std::string FileStart;
  
  const static int Precision;
  const static double CutoffCylinderLenth;
  const static double CutoffArrowLenth;
  const static double ArrowMaxConeRadiusFactor;
  const static double ArrowConeLengthFactor;
  
  std::string _tempFileName;
  std::ofstream _output;
  bool _currentlyInUnion;
  bool _currentlyInIntersectionUnion;
  double _aspectRatio;
};

#endif	/* POVRAYGENERATOR_H */

