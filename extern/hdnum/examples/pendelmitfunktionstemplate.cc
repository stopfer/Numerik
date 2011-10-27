// pendelmitfunktionstemplate.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen

template<typename Number>
void simuliere_pendel (Number l, Number phi, Number u)
{
  Number dt(1E-4);
  Number T(30.0);
  Number t(0.0);
  Number g(9.81/l);

  std::cout << t << " " << phi << std::endl;
  while (t<T)
  {
    t = t + dt;
    Number phialt(phi),ualt(u);
    phi = phialt + dt*ualt;
    u = ualt - dt*g*sin(phialt);
    std::cout << t << " " << phi << std::endl;
  }
}

int main () 
{
  float l1(1.34);  // Pendellänge in Meter
  float phi1(3.0); // Anfangsamplitude in Bogenmaß
  float u1(0.0);   // Anfangsgeschwindigkeit
  simuliere_pendel(l1,phi1,u1);

  double l2(1.34);  // Pendellänge in Meter
  double phi2(3.0); // Anfangsamplitude in Bogenmaß
  double u2(0.0);   // Anfangsgeschwindigkeit
  simuliere_pendel(l2,phi2,u2);
}
