// pendelnumerisch.cc
#include <iostream> // header für Ein-/Ausgabe
#include <cmath>    // mathematische Funktionen

int main ()
{
  double l(1.34);  // Pendellänge in Meter
  double phi(3.0); // Anfangsamplitude in Bogenmaß
  double u(0.0);   // Anfangsgeschwindigkeit
  double dt(1E-4); // Zeitschritt in Sekunden 
  double T(30.0);  // Ende in Sekunden
  double t(0.0);   // Anfangszeit

  std::cout << t << " " << phi << std::endl;
  while (t<T)
  {
    t = t + dt;        // inkrementiere Zeit
    double phialt(phi);// merke phi
    double ualt(u);    // merke u
    phi = phialt + dt*ualt;            // neues phi
    u = ualt - dt*(9.81/l)*sin(phialt);// neues u
    std::cout << t << " " << phi << std::endl;
  }
}
