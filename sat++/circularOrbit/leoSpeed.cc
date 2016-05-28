#include <iostream>
#include "constAstro.h"
#include "math.h"
#include "lcvProto.h"

// This is the main function, takes in a double precision altitude and returns the circular velocity
int main()
{
  using namespace constAstro;
  using namespace std;
  double SMA = getSMA();
  double vkps = getCircularSpeed(SMA);
  double vmph = kpsToMPH(vkps);
  cout << "\n";
  cout << "The LEO, circular speed of the satellite is " << vkps << " km/s" << endl;
  cout << "The LEO, circular speed of the satellite is " << vmph << " km/s" << endl;
  cout << "\n";

  return 0;
}
