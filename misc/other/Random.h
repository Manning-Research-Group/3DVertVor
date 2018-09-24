#ifndef RANDOM_H
#define RANDOM_H

#include "misc/geometry/Vector3D.h"

class Random {
public:
    // initialize random number generator...
    static unsigned int seed();
    static unsigned int seed(unsigned int s);
    
    // getter for seed
    static unsigned int getSeed() { return _seed; }
    
    // random number generators
    static int randomInteger();
    static double uniform(double maximum=1.0);
    static double exponential(double average=1.0);
    static double gamma(int shape, double average=1.0);
    static double normal(double standardDeviation=1.0, double average=0.0);
    static Vector3D gaussianVector3D(double standardDeviation=1.0);
    
private:
    static unsigned int _seed;
    static unsigned int getUrandom();
};

#endif
