#include <iostream>
#include <math.h>

double getSMA()
{
  std::cout << "\nEnter the altitude of the satellite in km: " << std::endl;
  double alt;              // Altitude in km
  double req = 6378.137;   // Earth radius in km
  std::cin >> alt;
  double SMA = alt + req;
  return SMA;
}

double getVelocity(double SMA, double mu)
{
  double v = sqrt(mu/SMA);
  return v;
}

int main()
{
  double mu = 3.986e5; // Gravitational parameter of Earth
  double SMA = getSMA();
  double v = getVelocity(SMA,mu);
  std::cout << "\n";
  std::cout << "The LEO, circular velocity of the satellite is: " << v << " km/s" << std::endl;
  std::cout << "\n";
  return 0;
}
