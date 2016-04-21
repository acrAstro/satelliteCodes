#include <iostream>
#include <math.h>

class TwoBodyOrbit
{
public:
  double J2;
  double initialKeplerElems[6][1];
  double safetyAltitude;
  double time;
};
