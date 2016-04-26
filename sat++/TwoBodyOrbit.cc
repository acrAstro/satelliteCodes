#include <iostream>
#include <math.h>

/* This is a class that does basic orbital analysis based on the Kepler elements
 * The primary functions are contained in the class, and the class is called
 * from a different file at compilation
*/
class TwoBodyOrbit
{
public:
  //double initialKeplerElems[6];
  //double period;
  //double meanMotion;
  double initTime;
  double finalTime;
  int numSteps;
  //double *time;
  void printTime()
  {
    using namespace std;
    cout  << "\n";
    cout << initTime << ":" << numSteps << ":" << finalTime << endl;
    cout << "\n";
  }
/*  double makeTimeVector()
  {
    time = new double[numSteps];
    return time;
  }*/
};

int main()
{
  double t0 = 0;
  double tf = 5400;
  double N = 200;

  TwoBodyOrbit tbp = {t0, tf, N};
  tbp.printTime();
  return 0;

}
