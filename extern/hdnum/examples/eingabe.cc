// eingabe.cc
#include <iostream> // header f�r Ein-/Ausgabe
#include <iomanip>  // f�r setprecision
#include <cmath>    // f�r sqrt
int main ()
{
  double x(0.0);
  std::cout << "Gebe eine Zahl ein: ";
  std::cin >> x;
  std::cout << "Wurzel(x)= " 
            << std::scientific << std::showpoint 
            << std::setprecision(15) 
            << sqrt(x) << std::endl;
}
