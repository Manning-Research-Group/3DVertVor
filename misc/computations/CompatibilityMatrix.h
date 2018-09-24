#ifndef COMPATIBILITYMATRIX_H
#define COMPATIBILITYMATRIX_H

#include <ostream>
#include <Eigen/Dense>

struct EigenDenseBackend {
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SelfAdjointEigenSolver<Matrix> Diagonalization;
  typedef Eigen::JacobiSVD<Matrix> Svd;
};

template <class Backend = EigenDenseBackend>
class CompatibilityMatrix {
public:
  typedef typename Backend::Vector Vector;
  typedef typename Backend::Matrix Matrix;
  typedef typename Backend::Svd Solver;
   
  CompatibilityMatrix(int numberOfConstraints, int numberOfDofs) : 
    _matrix(numberOfConstraints, numberOfDofs), 
    _solver(numberOfConstraints, numberOfDofs, Eigen::ComputeFullU | Eigen::ComputeFullV) {
    _matrix.fill(0.0);
  }
  int numberOfConstraints() const { return _matrix.rows(); }
  int numberOfDofs() const { return _matrix.cols(); }
  double operator()(int c, int dof) const { return _matrix(c,dof); }
  double &operator()(int c, int dof) { return _matrix(c,dof); }
  void addToElement(int c, int dof, double value) { _matrix(c,dof) += value; }
  friend std::ostream &operator<<(std::ostream &s, const CompatibilityMatrix &cm) { return s << cm._matrix; }
  void performSingularValueDecomposition() {
    _solver.compute(_matrix, Eigen::ComputeFullU | Eigen::ComputeFullV);
  }
  
  Matrix infinitesimalZeroModes(double cutoff=1e-10) {
    _solver.setThreshold(cutoff);
    const int Rank = _solver.rank();
    const Matrix &V = _solver.matrixV();
    return V.block(0,Rank,V.rows(),V.cols()-Rank);
  }
  Matrix infinitesimalSelfStresses(double cutoff=1e-10) {
    _solver.setThreshold(cutoff);
    const int Rank = _solver.rank();
    const Matrix &U = _solver.matrixU();
    return U.block(0,Rank,U.rows(),U.cols()-Rank);
  }
  const Vector &singularValues() const { return _solver.singularValues(); }
  const Matrix &u() const { return _solver.matrixU(); }
  const Matrix &v() const { return _solver.matrixV(); }
  const Matrix &matrix() const { return _matrix; }
  
private:
  Matrix _matrix;
  Solver _solver;
};

#endif /* COMPATIBILITYMATRIX_H */

