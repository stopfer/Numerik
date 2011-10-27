// pendelwhile.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen
int main ()
{
  double l(1.34);   // Pendellänge in Meter
  double phi0(0.2); // Amplitude im Bogenmaß
  double dt(0.05);  // Zeitschritt in Sekunden 
  double T(30.0);   // Ende in Sekunden
  double t(0.0);    // Anfangswert

  while ( t<=T )
  {
    std::cout << t << " " 
              << phi0*cos(sqrt(9.81/l)*t) 
              << std::endl;
    t = t + dt;
  }
}
