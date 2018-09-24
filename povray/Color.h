/* 
 * File:   Color.h
 * Author: arbeit
 *
 * Created on October 9, 2015, 9:56 PM
 */

#ifndef COLOR_H
#define	COLOR_H

#include <ostream>

class Color {
public:
  Color(double r, double g, double b) : _r(r), _g(g), _b(b) {}
  
  friend std::ostream& operator<<(std::ostream& os, const Color &c) { return os << c._r << ", " << c._g << ", " << c._b; }

private:
  double _r, _g, _b;
};

#endif	/* COLOR_H */

