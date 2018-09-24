#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Random.h"

unsigned int Random::_seed = Random::seed();

unsigned int Random::seed() {
  return seed(getUrandom());
}

unsigned int Random::seed(unsigned int s) {
  _seed = s;
  srand(s);
  printf("Using random seed 0x%08X\n", s);
  return s;
}

// by Doug; modified by Matthias:
unsigned int Random::getUrandom() {
  unsigned int randomSeed;
  const char fileName[] = "/dev/urandom";
  FILE *f = fopen(fileName, "rb");
  if(!f){
    printf("Failed to open random device %s.", fileName);
    return 0;
  }
  fread((char*)&randomSeed, sizeof(unsigned int), 1, f);
  fclose(f);
  return randomSeed;
}

// -------- RANDOM NUMBER GENERATORS ---------

int Random::randomInteger() {
  return rand();
}

double Random::uniform(double maximum) {
  return maximum*randomInteger()/((double)RAND_MAX+1);
}

double Random::exponential(double average) {
  return -log(1.0-uniform())*average;
}

double Random::gamma(int shape, double average) {
  double sum = 0.0;
  double scale = average/shape;
  if(shape<=0) {
    printf("Random::gamma: shape parameter must be >0, but is %d!", shape);
  }
  for(int i=0; i<shape; ++i) {
    sum += exponential(scale);
  }
  return sum;
}

// Box Muller method
double Random::normal(double standardDeviation, double average) {
  return average + standardDeviation*sqrt(exponential(2.0))*cos(uniform(2.0*M_PI));
}

Vector3D Random::gaussianVector3D(double standardDeviation) {
  return Vector3D(Random::normal(standardDeviation), Random::normal(standardDeviation), Random::normal(standardDeviation));
}
