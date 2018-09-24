#ifndef SPHERICAL_HARMONICS_FOR_VECTOR_H
#define SPHERICAL_HARMONICS_FOR_VECTOR_H

#include <cmath>
#include <complex>

#include "Vector3D.h"

/**
* Adopted from http://berenger.eu/blog/c-legendre-polynomial-by-recurrence-programming/  ( Berenger (contact at berenger dot eu) )
* This is the source code to construct the legendre polynome in C
* This is fast but you can improve the code by using pointer instead of
* accessing using index on the array and to compute (2*l-1) with a recurrence.
* Ref: Fast and accurate determination of the Wigner rotation matrices in FMM
*/
 
template <int MaxL> class SphericalHarmonicsForVector {
public:
  SphericalHarmonicsForVector(Vector3D v) {
    _psi = atan2(v.y(), v.x());
    const double z = v.z()/v.norm();
    const double sphericalHarmonicFactor = 0.5/sqrt(M_PI);
    const double sqrt1MinusZSq = sqrt(1.0-z*z);
 
    /**
    * Compute the legendre polynome by recurrence.
    * If needed you can use:
    * P_l,m (−x) = (−1)^m+l P_l,m(x)
    * Also:
    * P_l,-m (x) = (-1)^m * (l-m)!/(l+m)! * P_l,m(x)
    */    
    
    // Init legendre
    _legendre(0,0) = 1.0;
    _valuesWithoutPsiPart(0,0) = sphericalHarmonicFactor*_legendre(0,0);
    // Easy values
    const double factorSqrt3 = sphericalHarmonicFactor*sqrt(3);
    _legendre(1,0) = z;
    _legendre(1,1) = -sqrt1MinusZSq;
    _valuesWithoutPsiPart(1,0) = factorSqrt3*_legendre(1,0);
    _valuesWithoutPsiPart(1,1) = factorSqrt3*_legendre(1,1)/sqrt(2);
    _valuesWithoutPsiPart(1,-1) = -_valuesWithoutPsiPart(1,1);
 
    for(int l = 2; l <= MaxL; ++l) {
      double twoLMinus1 = 2*l-1;
      double factorSqrtTwoLPlus1 = sphericalHarmonicFactor*sqrt(2*l+1);
      _legendre(l,0) = (twoLMinus1*z*_legendre(l-1,0) - double(l-1)*_legendre(l-2,0))/double(l);
      _valuesWithoutPsiPart(l,0) = factorSqrtTwoLPlus1*_legendre(l,0);

      double facultyQuotient = 1.0;
      for(int m = 1; m < l-1; ++m) {
          // P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*P_l-2,m / (l-k)
          _legendre(l,m) = (twoLMinus1*z*_legendre(l-1,m) - double(l+m-1)*_legendre(l-2,m))/double(l-m);
          facultyQuotient /= double(l-m+1)*double(l+m);
          double sign = (m%2)?-1:1;
          _valuesWithoutPsiPart(l,m) = factorSqrtTwoLPlus1*_legendre(l,m)*sqrt(facultyQuotient);
          _valuesWithoutPsiPart(l,-m) = sign*factorSqrtTwoLPlus1*_legendre(l,m)*sqrt(facultyQuotient);
      }

      // P_l,l-1 = (2l-1)*x*P_l-1,l-1
      _legendre(l,l-1) = twoLMinus1*z*_legendre(l-1,l-1);
      double sign = (l%2)?-1:1;
      facultyQuotient /= 2*twoLMinus1;
      _valuesWithoutPsiPart(l,l-1) = factorSqrtTwoLPlus1*_legendre(l,l-1)*sqrt(facultyQuotient);
      _valuesWithoutPsiPart(l,-(l-1)) = -sign*factorSqrtTwoLPlus1*_legendre(l,l-1)*sqrt(facultyQuotient);

      // P_l,l = (2l-1)*factor*P_l-1,l-1
      _legendre(l,l) = -twoLMinus1*sqrt1MinusZSq*_legendre(l-1,l-1);
      facultyQuotient /= 2*l;
      _valuesWithoutPsiPart(l,l) = factorSqrtTwoLPlus1*_legendre(l,l)*sqrt(facultyQuotient);
      _valuesWithoutPsiPart(l,-l) = sign*factorSqrtTwoLPlus1*_legendre(l,l)*sqrt(facultyQuotient);
    }    
  }

  std::complex<double> value(int l, int m) const { return _valuesWithoutPsiPart(l,m)*std::complex<double>(cos(m*_psi),sin(m*_psi)); }
private:
  double &_legendre(int l, int m) { return _legendreArray[l][m]; }
  double _valuesWithoutPsiPart(int l, int m) const { return _valuesWithoutPsiPartArray[l][l+m]; }
  double &_valuesWithoutPsiPart(int l, int m) { return _valuesWithoutPsiPartArray[l][l+m]; }
  double _legendreArray[MaxL+1][MaxL+1];
  double _valuesWithoutPsiPartArray[MaxL+1][2*MaxL+1];
  double _psi;
};

//
//void computeLegendre(FReal legendre[], const FReal x, const int P){
//    //This factor is reuse −sqrt(1 − x^2)
//    const FReal factor = -sqrt(1.0-pow(x,2));
// 
//    // Init legendre
//    legendre[atLm(0,0)] = 1.0;        // P_0,0(x) = 1
//    // Easy values
//    legendre[atLm(1,0)] = x;      // P_1,0(x) = x
//    legendre[atLm(1,1)] = factor;     // P_1,1(x) = −sqrt(1 − x^2)
// 
//    for(int l = 2; l <= P ; ++l ){
//        for( int m = 0; m < l - 1 ; ++m ){
//            // P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k)
//            legendre[atLm(l,m)] = (FReal(2*l-1) * x * legendre[atLm(l-1,m)] - FReal( l + m - 1 ) * legendre[atLm(l-2,m)] )
//                    / FReal( l - m );
//        }
//        // P_l,l-1 = (2l-1)*x*P_l-1,l-1
//        legendre[atLm(l,l-1)] = FReal(2*l-1) * x * legendre[atLm(l-1,l-1)];
//        // P_l,l = (2l-1)*factor*P_l-1,l-1
//        legendre[atLm(l,l)] = FReal(2*l-1) * factor * legendre[atLm(l-1,l-1)];
//    }
//}
 
///**
//* To test the legendre function
//*/
//int main(){
//    static const int P = 2;
//    FReal legendre[((P+2)*(P+1))/2];
// 
//    computeLegendre(legendre, 0.5, P);
// 
//    for(int l = 0 ; l <= P ; ++l){
//        for(int m = 0 ; m <= l ; ++m){
//            printf("P{%d,%d} = %lf t", l, m, legendre[atLm(l,m)]);
//        }
//        printf("n");
//    }
// 
//    return 0;
//}

#endif /* LEGENDEPOLYNOMIALS_H */

