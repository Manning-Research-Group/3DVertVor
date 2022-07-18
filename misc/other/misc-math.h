#ifndef MISC_MATH_H
#define	MISC_MATH_H

#include <math.h>

inline double cubicRoot(const double x) {
  return exp(log(x)/3.0);
}

inline int modulo(int i, int n) {
  return (i % n + n) % n;
}

#endif	/* MISC_MATH_H */

