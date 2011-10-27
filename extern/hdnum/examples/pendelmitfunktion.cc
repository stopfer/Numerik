// pendelmitfunktion.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen

void simuliere_pendel (double l, double phi, double u) 
{
  double dt   = 1E-4;
  double T    = 30.0;
  double t    = 0.0;

  std::cout << t << " " << phi << std::endl;
  while (t<T) 
  {
    t = t + dt;
    double phialt(phi),ualt(u);
    phi = phialt + dt*ualt;
    u = ualt - dt*(9.81/l)*sin(phialt);
    std::cout << t << " " << phi << std::endl;
  }
}

int main () 
{
  double l(1.34); // Pendellänge in Meter
  double phi(3.0);// Anfangsamplitude in Bogenmaß
  double u(0.0);  // Anfangsgeschwindigkeit
  simuliere_pendel(l,phi,u);
}
