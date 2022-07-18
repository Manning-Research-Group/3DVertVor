#ifndef DYNAMICALMATRIX_H
#define DYNAMICALMATRIX_H

#include <vector>

#include <Eigen/Dense>

class Hessian : public Eigen::MatrixXd {
public:
  Hessian(int numberOfDofs) : Eigen::MatrixXd(numberOfDofs, numberOfDofs) { fill(0.0); }
  int numberOfDofs() const { return rows(); }
};

#endif /* DYNAMICALMATRIX_H */

