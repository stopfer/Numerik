// pendelmittimer.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen
#include "hdnum.hh" // Zeitmessung
int main ()
{
  double l(1.34);   // Pendellänge in Meter
  double phi0(0.2); // Amplitude im Bogenmaß
  double dt(0.05);  // Zeitschritt in Sekunden 
  double T(30.0);   // Ende in Sekunden
  hdnum::Timer zeit,zeitIter;
  for (double t=0.0; t<=T; t=t+dt)
  {
    zeitIter.reset();
    std::cout << t << " " 
              << phi0*cos(sqrt(9.81/l)*t) 
              << std::endl;
    std::cout << "Durchgang " << int(t/dt) << " brauchte "
              << std::scientific << zeitIter.elapsed()
              << " Sekunden" << std::endl;
  }
  std::cout << "Die Ausgabe aller Werte brauchte " << zeit.elapsed()
            << " Sekunden" << std::endl;
}
