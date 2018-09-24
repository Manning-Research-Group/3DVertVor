#ifndef DYNAMICALMATRIXOLD_H
#define DYNAMICALMATRIXOLD_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class DynamicalMatrixOld {
public:
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SelfAdjointEigenSolver<Matrix> Solver;
  
  DynamicalMatrixOld(int numberOfPositionScalars, int numberOfControlParameters)
      : _numberOfPositionScalars(numberOfPositionScalars), _numberOfControlParameters(numberOfControlParameters),
      _mixedDerivativesEnergyPositionPosition(_numberOfPositionScalars, _numberOfPositionScalars), _eigensolver(NULL) {
    reset();
  }
  virtual ~DynamicalMatrixOld() {
    cleanup();
  }
  void reset() {
    cleanup();
    _mixedDerivativesEnergyPositionPosition.fill(0.0);
    for(int i=0; i<_numberOfControlParameters; ++i) {
      _mixedDerivativesEnergyControlParameterPosition.push_back(Vector(_numberOfPositionScalars));
      _mixedDerivativesEnergyControlParameterPosition.back().fill(0.0);
      _2ndDerivativeControlParameter.push_back(0.0);
    }
  }
  void addDerivativePositionPosition(int p1, int p2, double derivative) {
    _mixedDerivativesEnergyPositionPosition(p1,p2) += derivative;
  }
  void addDerivativeControlParameterPosition(int c, int p, double derivative) {
    _mixedDerivativesEnergyControlParameterPosition[c](p) += derivative;
  }
  void add2ndDerivativeControlParameter(int c, double derivative) {
    _2ndDerivativeControlParameter[c] += derivative;
  }
  DynamicalMatrixOld &operator*=(double factor) {
    _mixedDerivativesEnergyPositionPosition *= factor;
    for(int i=0; i<_numberOfControlParameters; ++i) {
      _mixedDerivativesEnergyControlParameterPosition[i] *= factor;
      _2ndDerivativeControlParameter[i] *= factor;
    }
    return *this;
  }
  DynamicalMatrixOld &operator+=(const DynamicalMatrixOld &other) {
    assert(_numberOfPositionScalars==other._numberOfPositionScalars);
    assert(_numberOfControlParameters==other._numberOfControlParameters);
    _mixedDerivativesEnergyPositionPosition += other._mixedDerivativesEnergyPositionPosition;
    for(int i=0; i<_numberOfControlParameters; ++i) {
      _mixedDerivativesEnergyControlParameterPosition[i] += other._mixedDerivativesEnergyControlParameterPosition[i];
      _2ndDerivativeControlParameter[i] += other._2ndDerivativeControlParameter[i];
    }
    return *this;
  }
  bool solve() {
    _eigensolver = new Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(_mixedDerivativesEnergyPositionPosition);
    return _eigensolver->info() == Eigen::Success;
  }
  const Solver::RealVectorType& eigenValues() const {
    return _eigensolver->eigenvalues();
  }
  const Solver::EigenvectorsType& eigenVectors() const {
    return _eigensolver->eigenvectors();
  }
  void printEigenValuesTo(std::ostream &os, bool onlyValues=false) const {
    for(int i=0; i<_eigensolver->eigenvalues().rows(); ++i) {
      if(onlyValues) {
        os << _eigensolver->eigenvalues()(i) << std::endl;
      } else {
        os << i << " .. " << _eigensolver->eigenvalues()(i) << std::endl;
      }
    }
  }
  void printEigenValues() const { printEigenValuesTo(std::cout); }
  double modulusWithRespectToControlParameter(int c, double eigenValueCutoff) { // normal one
    double modulus = _2ndDerivativeControlParameter[c];
    for(int i=0; i<eigenValues().rows(); ++i) {
      double omegaSq = eigenValues()(i);
      if(omegaSq>=eigenValueCutoff) {
        double projection = eigenVectors().col(i).transpose() * _mixedDerivativesEnergyControlParameterPosition[c];
        modulus -= projection*projection/omegaSq;
      }
    }
    return modulus;
  }
  double scalarProductEigenvectorMixedControlParameterDerivative(int cp, int i) const { return eigenVectors().col(i).transpose() * _mixedDerivativesEnergyControlParameterPosition[cp]; }
  double operator()(int i1, int i2) const { return _mixedDerivativesEnergyPositionPosition(i1,i2); }
  double mixedDerivative(int cp, int i) const { return _mixedDerivativesEnergyControlParameterPosition[cp](i); }
  double controlParameterDerivative(int cp) const { return _2ndDerivativeControlParameter[cp]; }
  const Vector &mixedDerivatives(int cp) const { return _mixedDerivativesEnergyControlParameterPosition[cp]; }
  const Matrix &positionHessian() const { return _mixedDerivativesEnergyPositionPosition; }
  
private:
  void cleanup() {
    if(_eigensolver) {
      delete _eigensolver;
      _eigensolver = NULL;
    }
    _mixedDerivativesEnergyControlParameterPosition.clear();
    _2ndDerivativeControlParameter.clear();
  }
  
  int _numberOfPositionScalars;
  int _numberOfControlParameters;
  Matrix _mixedDerivativesEnergyPositionPosition;
  std::vector<Vector> _mixedDerivativesEnergyControlParameterPosition;
  std::vector<double> _2ndDerivativeControlParameter;
  Solver *_eigensolver;
};

#endif /* DYNAMICALMATRIX_H */

