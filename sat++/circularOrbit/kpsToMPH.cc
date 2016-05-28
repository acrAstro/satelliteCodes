#include "constAstro.h"

// This function converts the circular velocity from km/s to mph
double kpsToMPH(double kps)
{
  using namespace constAstro;
  double mph = kps*cf;
  return mph;
}
