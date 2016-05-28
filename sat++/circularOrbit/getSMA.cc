#include <iostream>
#include "constAstro.h"
#include "math.h"

/*
This function takes the altitude of a satellite in a circular orbit and
returns the semi-major axis
*/
double getSMA()
{
  using namespace std;
  using namespace constAstro;
  cout << "\nEnter the altitude of the satellite in km: " << endl;
  double alt;
  cin >> alt;
  double SMA = alt + req;
  return SMA;
}
