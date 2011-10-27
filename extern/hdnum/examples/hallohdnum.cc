// hallohdnum.cc 
#include <iostream>    // notwendig zur Ausgabe
#include <vector>
#include "hdnum.hh"    // hdnum header

int main ()
{
  hdnum::Vector<float> a(10,3.14); // Feld mit 10 init. Elementen
  a[3] = 1.0;              // Zugriff auf Element 3
}
