#include "Matrix3x3x3.h"

const Matrix3x3x3 Matrix3x3x3::Zero = Matrix3x3x3(Matrix3x3::Zero, Matrix3x3::Zero, Matrix3x3::Zero);

const Matrix3x3x3 Matrix3x3x3::Epsilon = Matrix3x3x3(Matrix3x3(Vector3D(0,0,0),Vector3D(0,0,-1),Vector3D(0,1,0)),
                                                     Matrix3x3(Vector3D(0,0,1),Vector3D(0,0,0),Vector3D(-1,0,0)),
                                                     Matrix3x3(Vector3D(0,-1,0),Vector3D(1,0,0),Vector3D(0,0,0)));
