/* 
 * File:   PovrayGenerator.cpp
 * Author: arbeit
 * 
 * Created on October 9, 2015, 9:35 PM
 */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <unordered_map>

#include "PovrayGenerator.h"
#include "misc/other/fileIo.h"

const std::string PovrayGenerator::TempFolderName("temp");
const std::string PovrayGenerator::TempFileNameMask("temp/povXXXXXX");
const std::string PovrayGenerator::LogFileName("temp/povray.log");
/*const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
background{rgb <1.00,1.00,1.00>}\n\
light_source{<8,-20,30> color rgb <0.77,0.75,0.75>}\n\
light_source{<-8,20,30> color rgb <0.77,0.75,0.75>}\n\
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}\n\
");
const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
background{rgb <1.00,1.00,1.00>}\n\
light_source{<8,-20,30> color rgb <0.77,0.75,0.75>}\n\
light_source{<-8,20,30> color rgb <0.77,0.75,0.75>}\n\
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}\n\
");
const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
background{rgb <0.85,0.90,1.00>}\n\
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}\n\
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}\n\
");
const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
background{rgb <0.92,1.00,0.90>}\n\
light_source{<200,0,0> color rgb <0.75,0.70,0.70>}\n\
light_source{<0,-200,0> color rgb <0.30,0.30,0.30>}\n\
light_source{<0,0,200> color rgb <0.95,0.90,0.90>}\n\
");
const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
background{rgb <1.00,1.00,1.00>}\n\
light_source{<200,0,0> color rgb <0.75,0.70,0.70>}\n\
light_source{<0,-200,0> color rgb <0.30,0.30,0.30>}\n\
light_source{<0,0,200> color rgb <0.95,0.90,0.90>}\n\
");*/
const std::string PovrayGenerator::FileStart("\
#version 3.6;\n\
");

const int PovrayGenerator::Precision = 8;
const double PovrayGenerator::CutoffCylinderLenth = 1e-6;
const double PovrayGenerator::CutoffArrowLenth = 1e-6;
const double PovrayGenerator::ArrowMaxConeRadiusFactor = 2.4;
const double PovrayGenerator::ArrowConeLengthFactor = 5.0;


PovrayGenerator::PovrayGenerator() : _currentlyInUnion(false), _currentlyInIntersectionUnion(false), _aspectRatio(0.75) {
}

PovrayGenerator::~PovrayGenerator() {
}

void PovrayGenerator::start(const Color &background) {
  recursivelyCreateFolder(TempFolderName);
  char temp[256];
  strcpy(temp, TempFileNameMask.c_str());
  int fd = mkstemp(temp);
  if(fd!=-1) {
    close(fd);
    _tempFileName = temp;
//    std::cout << "Preparing Povray file " << _tempFileName << "..." << std::endl;
    _output.open(_tempFileName.c_str());
    _output.precision(Precision);
    _output << FileStart;
    _output << "background{rgb <" << background << ">}\n";
  } else {
    std::cout << "Could not create temporary file for povray source!" << std::endl;
    exit(1);
  }
}

void PovrayGenerator::create(std::string filename, int width) {
  // end last union if that has been forgotten
  endUnion();
  // close file
  _output.close();
  
  // run povray
  std::cout << "Rendering to " << filename << "..." << std::endl;
  std::stringstream ss;
  ss << "povray -D +W" << width << " +H" << round(width/_aspectRatio) << " +A0.3 +O" << filename << " +I" << _tempFileName << " 2>" << LogFileName;
  std::cout << ss.str() << std::endl;
  if(0==std::system(ss.str().c_str())) {
    std::cout << "Done." << std::endl;  
    if(remove(_tempFileName.c_str())) {
      std::cout << "Error removing " << _tempFileName  << "!" << std::endl;
    }
  } else {
    std::cout << "Error rendering!" << std::endl;
  }
}

void PovrayGenerator::addLight(const Vector3D &position, const Color &color) {
  _output << "light_source{<" << position << "> color rgb <" << color << ">}\n";
}


void PovrayGenerator::setCamera(const Vector3D &position, const Vector3D &looksAt, const double aspectRatio) {
  _aspectRatio = aspectRatio;
  _output << "camera {\n";
  _output << "  location <" << position << ">\n";
  _output << "  up <0,1.0,0.0>\n";
  _output << "  right <" << -_aspectRatio << ",0,0>\n";
  _output << "  sky z\n";
  _output << "  look_at <" << looksAt << ">\n";
  _output << "}\n";
}

void PovrayGenerator::setVariable(const std::string &name, double value) {
  _output << "#declare " << name << " = " << value << ";\n";
}


bool PovrayGenerator::startUnion() {
  if(_currentlyInUnion) {
    return false;
  }
  
  _output << "union{\n";
  _currentlyInUnion = true;
  return true;
}

bool PovrayGenerator::startUnionWithBoxIntersection(const Vector3D &Corner1, const Vector3D &Corner2) {
  if(_currentlyInUnion) {
    return false;
  }
  
  _output << "union{intersection{\n";
  _output << "  box{<" << Corner1 << ">, <" << Corner2 << ">}\n";
  _output << "  union{\n";
  _currentlyInUnion = true;
  _currentlyInIntersectionUnion = true;
  return true;
}

void PovrayGenerator::endUnion() {
  if(_currentlyInUnion) {
    if(_currentlyInIntersectionUnion) {
      _output << "}}}\n";
    } else {
      _output << "}\n";
    }
    _currentlyInUnion = false;
    _currentlyInIntersectionUnion = false;
  }
}

void PovrayGenerator::endUnionWithColor(const Color &color) {
  if(_currentlyInUnion) {
    if(_currentlyInIntersectionUnion) {
      _output << "  }}\n";
    }    
    _output << "  pigment{rgb <" << color << ">} finish{specular 0.5 ambient 0.42}\n";
    _output << "}\n";
    _currentlyInUnion = false;
    _currentlyInIntersectionUnion = false;
  }
}

void PovrayGenerator::endUnionWithColorF(const Color &color, const double filterTransparency) {
  if(_currentlyInUnion) {
    if(_currentlyInIntersectionUnion) {
      _output << "  }}\n";
    }    
    _output << "  pigment{rgbf <" << color << "," << filterTransparency << ">} finish{specular 0.5 ambient 0.42}\n";
    _output << "}\n";
    _currentlyInUnion = false;
    _currentlyInIntersectionUnion = false;
  }
}

void PovrayGenerator::endUnionWithColorT(const Color &color, const double transmTransparency) {
  if(_currentlyInUnion) {
    if(_currentlyInIntersectionUnion) {
      _output << "  }}\n";
    }    
    _output << "  pigment{rgbt <" << color << "," << transmTransparency << ">} finish{specular 0.5 ambient 0.42}\n";
    _output << "}\n";
    _currentlyInUnion = false;
    _currentlyInIntersectionUnion = false;
  }
}


void PovrayGenerator::drawCoordinateSystem() {
  startUnion();
  drawArrow(Vector3D(0,0,0), Vector3D(1,0,0), 0.05);
  endUnionWithColor(Color(1.0,0.0,0.0));

  startUnion();
  drawArrow(Vector3D(0,0,0), Vector3D(0,1,0), 0.05);
  endUnionWithColor(Color(0.0,1.0,0.0));

  startUnion();
  drawArrow(Vector3D(0,0,0), Vector3D(0,0,1), 0.05);
  endUnionWithColor(Color(0.0,0.0,1.0));
}

void PovrayGenerator::drawSphere(const Vector3D &c, const double radius) {
  _output << "  sphere{<" << c << ">," << radius << "}\n";
}

void PovrayGenerator::drawCylinder(const Vector3D &p1, const Vector3D &p2, const double radius) {
  if((p2-p1).norm()>CutoffCylinderLenth) {
    _output << "  cylinder{<" << p1 << ">,<" << p2 << ">," << radius << "}\n";
  }
}

void PovrayGenerator::drawArrow(const Vector3D &p1, const Vector3D &p2, const double radius) {
  double totalLength = (p2-p1).norm();
  if(totalLength>CutoffArrowLenth) {
    double maxConeRadius = ArrowMaxConeRadiusFactor*radius;
    double coneLength = ArrowConeLengthFactor*radius;
    Vector3D unitVector((p2-p1)/totalLength);
    if(totalLength>coneLength+CutoffCylinderLenth) {
      _output << "  cylinder{<" << p1 << ">,<" << p1 + (totalLength-coneLength)*unitVector << ">," << radius << "}\n";
    }
    _output << "  cone{<" << p2-coneLength*unitVector << ">," << maxConeRadius << "\n";
    _output << "       <" << p2 << ">,0}\n";
  }
}
