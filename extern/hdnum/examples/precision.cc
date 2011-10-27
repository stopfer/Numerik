// precision.cc (Dateiname als Kommentar)
#include <iostream>  // notwendig zur Ausgabe
#include <iomanip>   // für setprecision etc.
#include "hdnum.hh"

int main ()
{
  float eps_float;
  int i_float = hdnum::precision(eps_float);
  std::cout << "float: " << "2^-" << i_float << std::endl;
  double eps_double;
  int i_double = hdnum::precision(eps_double);
  std::cout << "double: " << "2^-" << i_double << std::endl;
}
