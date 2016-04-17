#include <iostream>
#include <math.h>

/* This function takes the altitude of a satellite in a circular orbit and
 * returns the semi-major axis */
double getSMA()
{
  std::cout << "\nEnter the altitude of the satellite in km: " << std::endl;
  double alt;              // Altitude in km
  double req = 6378.137;   // Earth radius in km
  std::cin >> alt;
  double SMA = alt + req;
  return SMA;
}

/* This function computes the velocity of the satellite in km/s from the
 * semi-major axis */
double getVelocity(double SMA, double mu)
{
  double v = sqrt(mu/SMA);
  return v;
}

// This function converts the circular velocity from km/s to mph
double kpsToMPH(double kps)
{
  double cf = 2236.936; // 1 km/s is 2236.936 mph
  double mph = kps*cf;
  return mph;
}

// Main function!!!
int main()
{
  double mu = 3.986e5; // Gravitational parameter of Earth
  double SMA = getSMA();
  double vkps = getVelocity(SMA,mu);
  double vmph = kpsToMPH(vkps);
  std::cout << "\n";
  std::cout << "The LEO, circular velocity of the satellite is: " << vkps << " km/s" << std::endl;
  std::cout << "The LEO, circular velocity of the satellite is: " << vmph << " mph"  << std::endl;
  std::cout << "\n";

  return 0;
}
