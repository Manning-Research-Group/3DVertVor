#ifndef LMDERIVATIVEBISECTION_H
#define LMDERIVATIVEBISECTION_H

#include <iostream>
#include <functional>
#include <math.h>
#include <sstream>
#include <exception>

class LmDerivativeNewton {
public:
  LmDerivativeNewton() : _c2(1e-4), _minMaxAbsSlope(1e-10), _maxInitialAlphaFactor(2), _factorSearchForMaxAlpha(1.5), _maxIterations(100), _debug(false), _iterations(0) {}
  void setC2(double c2) { _c2 = c2; }
  void setMinMaxAbsSlope(double minMaxAbsSlope) { _minMaxAbsSlope = minMaxAbsSlope; }
  void debug() { _debug=true; }

  long iterations() const { return _iterations; }
  
  double minimize(double lastAlphaTimesLastSlope, double initialSlope, const std::function<double(double)> &functionValue, const std::function<double(double)> &gradient, const std::function<double(double,double&)> &functionValueAndGradient) {
    // abort criterion
    double mas = _c2*fabs(initialSlope);
    const double MaxAbsSlope = (mas>_minMaxAbsSlope)?mas:_minMaxAbsSlope;
    const double InitialAlpha = lastAlphaTimesLastSlope/fabs(initialSlope);
    const double MaxSecondAlpha = _maxInitialAlphaFactor*InitialAlpha;

    
    struct Point {
      Point(double a, const std::function<double(double)> &gradient) : alpha(a), slope(gradient(alpha)) {/*debugPrint("constructor");*/}
      Point(double a, double s) : alpha(a), slope(s) {/*debugPrint("constructor");*/}
      void setAlpha(double a, const std::function<double(double)> &gradient) { alpha=a; slope = gradient(alpha); /*debugPrint("setAlpha");*/}
      double alpha, slope;
      bool check(const double MaxAbsSlope) const { return fabs(slope)<MaxAbsSlope; }
      bool operator<(const Point &o) const { return alpha < o.alpha; }
      bool operator>(const Point &o) const { return alpha > o.alpha; }
      Point operator-(const Point &o) { Point tmp(*this); tmp.alpha-=o.alpha; tmp.slope-=o.slope; return tmp; }
      void swap(Point &o) { Point temp = *this; *this = o; o = temp; }
//      void debugPrint(const std::string &context="") const { 
//        std::cout << context << ":  position: " << alpha << " -> slope: " << slope << "\n"; 
//      }
    };

    // first two points
    Point minimum(0, initialSlope);
    Point middle(InitialAlpha, gradient);

    // maybe we are lucky
    if(middle.check(MaxAbsSlope)) return middle.alpha;

    // linear inter-/extrapolation of slope to get alpha max
    double maxAlpha = fabs(middle.alpha*minimum.slope/(middle.slope-minimum.slope));
    if(maxAlpha>MaxSecondAlpha) maxAlpha = MaxSecondAlpha;
    Point maximum(maxAlpha, gradient);
//    Point maximum( fabs(middle.alpha*minimum.slope/(middle.slope-minimum.slope)), gradient);

    // maybe we are lucky
    if(maximum.check(MaxAbsSlope)) return maximum.alpha;
    
    // check order:
    if(middle > maximum) middle.swap(maximum);
    
    // make sure we have a t least one point with slope>0 (middle or maximum))
    if(middle.slope<0) {
      while(maximum.slope<0) {
        double maxAlpha = maximum.alpha;
        minimum = middle;
        middle = maximum;
        maximum.setAlpha(maxAlpha*_factorSearchForMaxAlpha, gradient);
        if(maximum.check(MaxAbsSlope)) return maximum.alpha;
      };
    }
    
    _iterations = 0;
    while(true) {
      if(_iterations>_maxIterations) {
        std::cout << "LmDerivativeNewton::minimize:  maximum number of iterations!\n";
        return middle.alpha;
      }
      
//      std::cout << "Bla.\n";
      // quadratic interpolation of slope;  parabola relative to minimum.alpha:
      Point deltaMid = middle-minimum;
      double deltaAlphaSqMid = deltaMid.alpha*deltaMid.alpha;
      Point deltaMax = maximum-minimum;
      double deltaAlphaSqMax = deltaMax.alpha*deltaMax.alpha;
      double denominator = deltaMid.alpha*deltaAlphaSqMax - deltaMax.alpha*deltaAlphaSqMid;
      if(denominator==0) {
//        std::cout << "LmDerivativeNewton::minimize:  zero denominator!  Slope: " << middle.slope << "\n";
        return middle.alpha;
      }
      double c = (deltaMax.slope*deltaMid.alpha - deltaMid.slope*deltaMax.alpha)/denominator;
      double b = (deltaMid.slope*deltaAlphaSqMax - deltaMax.slope*deltaAlphaSqMid)/denominator;
      double a = minimum.slope;
      
      // compute zero:
      double offset = 0.5*b/c;
      double det = offset*offset - a/c;
      if(det<-1e-14*offset*offset) {
//        std::cout << "LmDerivativeNewton::minimize:  Got negative determinant: " << det << ", quotient: " << det/(offset*offset) << "\n";
      }
      double sqrtDet = (det<=0)?0:sqrt(det);
//      std::cout << "c = " << c << "\n";
      double newAlpha = (c>0) ? -offset + sqrtDet : -offset - sqrtDet;
//      std::cout << "newAlpha: " << newAlpha << "\n";
//      std::cout << "check: " << a + b*newAlpha + c*newAlpha*newAlpha << "\n";
      Point subMiddle(minimum.alpha + newAlpha, gradient);
      if(subMiddle.check(MaxAbsSlope)) {
//        std::cout << "lm slope: " << subMiddle.slope << ", initial:  " << initialSlope << ", it:  " << _iterations << "\n";
        return subMiddle.alpha;
      }

      if(_debug) {
        std::cout << "checkMid: " << a + b*middle.alpha + c*middle.alpha*middle.alpha << "\n";
        std::cout << "checkMax: " << a + b*maximum.alpha + c*maximum.alpha*maximum.alpha << "\n";
        std::cout << "cutoff slope:  " << MaxAbsSlope << "\n";
        std::cout << "c:  " << c << "\n";
        std::cout << "minimum.alpha*c:  " << minimum.alpha*c << "\n";
        std::cout << "b:  " << b << "\n";
        std::cout << "minimum.alpha*(b + minimum.alpha*c):  " << minimum.alpha*(b + minimum.alpha*c) << "\n";
        std::cout << "minimum.slope:  " << minimum.slope << "\n";
        std::cout << "a:  " << a << "\n";
        std::cout << "offset:  " << offset << "\n";
        std::cout << "offset*offset:  " << offset*offset << "\n";
        std::cout << "a/c:  " << a/c << "\n";
        std::cout << "sqrtDet:  " << sqrtDet << "\n";
        std::cout << "check: " << a + b*newAlpha + c*newAlpha*newAlpha << "\n";
        std::cout << "min abs:  " << minimum.alpha << ", slope: " << minimum.slope << "\n";
        std::cout << "smid: " << subMiddle.alpha-minimum.alpha << ", slope: " << subMiddle.slope  << ", rel. slope: " << subMiddle.slope-minimum.slope  << "\n";
        std::cout << "mid:  " << middle.alpha-minimum.alpha << ", slope: " << middle.slope << ", rel. slope: " << middle.slope-minimum.slope << "\n";
        std::cout << "max:  " << maximum.alpha-minimum.alpha << ", slope: " << maximum.slope << ", rel. slope: " << maximum.slope-minimum.slope << "\n";
      }
      
      // sort alphas:
      if(subMiddle > middle) subMiddle.swap(middle);

      // decide which one to kick out (either minimum or maximum):
      if(subMiddle.slope>0) {
        // kick out maximum, because we want minimum.slope always <0
        maximum = middle;
        middle = subMiddle;
      } else if(middle.slope<0) {
        // kick out minimum, because we want one slope >0, which is maximum
        minimum = subMiddle;
      } else {
        // two above and two below -> compare slopes
        if(fabs(minimum.slope)<fabs(maximum.slope)) {
          // minimum is closer to zero, so kick out maximum
          maximum = middle;
          middle = subMiddle;
        } else {
          // maximum is closer to zero, so kick out minimum
          minimum = subMiddle;
        }
      }
      ++_iterations;
    }
    
    return middle.alpha;
  }
  
private:
  double _c2;  // second, strong wolfe condition
  double _minMaxAbsSlope;
  double _maxInitialAlphaFactor;
  double _factorSearchForMaxAlpha;
  long _maxIterations;
  bool _debug;
  long _iterations;
};

#endif /* LMDERIVATIVEBISECTION_H */

