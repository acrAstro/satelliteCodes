#include <iostream>
#include "math.h"
#include "constAstro.h"

/*
This function computes the velocity of the satellite in km/s from the semi-major axis
*/
double getCircularSpeed(double SMA)
{
  using namespace constAstro;
  double v = sqrt(mu/SMA);
  return v;
}
