#ifndef DYNAMICALMATRIXTEST_H
#define DYNAMICALMATRIXTEST_H

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class DynamicalMatrixTest {
public:
  typedef Eigen::VectorXd Vector;
  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SelfAdjointEigenSolver<Matrix> Solver;
  
  DynamicalMatrixTest(int numberOfPositionScalars, int numberOfControlParameters /* will be set to one - this is just a test */ )
      : _numberOfPositionScalars(numberOfPositionScalars), _numberOfControlParameters(1),
      _secondDerivatives(_numberOfPositionScalars+_numberOfControlParameters, _numberOfPositionScalars+_numberOfControlParameters), _eigensolver(NULL) {
    reset();
  }
  virtual ~DynamicalMatrixTest() {
    cleanup();
  }
  void reset() {
    cleanup();
    _secondDerivatives.fill(0.0);
  }
  void addDerivativePositionPosition(int p1, int p2, double derivative) {
    _secondDerivatives(p1,p2) += derivative;
  }
  void addDerivativeControlParameterPosition(int c, int p, double derivative) {
    if(c==0) {
      _secondDerivatives(p,_numberOfPositionScalars+c) += derivative;
      _secondDerivatives(_numberOfPositionScalars+c,p) = _secondDerivatives(p,_numberOfPositionScalars+c);
    }
  }
  void add2ndDerivativeControlParameter(int c, double derivative) {
    if(c==0) {
      _secondDerivatives(_numberOfPositionScalars+c,_numberOfPositionScalars+c) += derivative;
    }
  }
  bool solve() {
    _eigensolver = new Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(_secondDerivatives);
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
    double sum = 0.0;
    for(int i=0; i<eigenValues().rows(); ++i) {
      double omegaSq = eigenValues()(i);
//      if(omegaSq>=eigenValueCutoff) {
        double element = eigenVectors()(_numberOfPositionScalars+c, i);
        sum += element*element/omegaSq;
//      } else {
//        if()
//      }
    }
    return 1.0/sum;
  }
//  double scalarProductEigenvectorMixedControlParameterDerivative(int cp, int i) const { return eigenVectors().col(i).transpose() * _mixedDerivativesEnergyControlParameterPosition[cp]; }
//  double operator()(int i1, int i2) const { return _secondDerivatives(i1,i2); }
//  double mixedDerivative(int cp, int i) const { return _mixedDerivativesEnergyControlParameterPosition[cp](i); }
//  double controlParameterDerivative(int cp) const { return _2ndDerivativeControlParameter[cp]; }
//  Matrix positionHessian() const { return _secondDerivatives; }
  
private:
  void cleanup() {
    if(_eigensolver) {
      delete _eigensolver;
      _eigensolver = NULL;
    }
  }
  
  int _numberOfPositionScalars;
  int _numberOfControlParameters;
  Matrix _secondDerivatives;
  Solver *_eigensolver;
};

#endif /* DYNAMICALMATRIX_H */

