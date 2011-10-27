// vektoren.cc 
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

template<class T>
void product (hdnum::Vector<T> &x)
{
  for (int i=1; i<x.size(); i=i+1)
    x[i] = x[i]*x[i-1];
}

template<class T>
T sum (hdnum::Vector<T> x)
{
  T s(0.0);
  for (int i=0; i<x.size(); i=i+1)
    s = s + x[i];
  return s;
}

int main ()
{
  // Konstruktion
  hdnum::Vector<float> x(10);        // Vektor mit 10 Elementen
  hdnum::Vector<double> y(10,3.14);  // 10 Elemente initialisiert
  hdnum::Vector<float> a;            // ein leerer Vektor
  x.resize(117);              // vergrößern, Daten gelöscht!
  x.resize(23,2.71);          // verkleinern geht auch

  // Zugriff auf Vektorelemente
  for (int i=0; i<x.size(); i=i+1)
  x[i] = i;                 // Zugriff auf Elemente 

  // Kopie und Zuweisung
  hdnum::Vector<float> z(x);         // Kopie erstellen
  z[2] = 1.24;                // Wert verändern
  
  a = z;                      // hat Werte von z
  a[2] = -0.33;               
  a = 5.4;                    // Zuweisung an alle Elemente
  
  hdnum::Vector<float> w(x);        
  
  // Rechnen mit Vektoren
  w += z;                     // w = w+z
  w -= z;                     // w = w-z
  w *= 1.23;                  // skalare Multiplikation
  w /= 1.23;                  // skalare Division
  w.update(1.23,z);           // w = w + a*z
  x[0] = w*z;                 // skalare Multiplikation
  std::cout << x.two_norm() << std::endl; // euklidische Norm

  // Ausgabe
  std::cout << w << std::endl;// schöne Ausgabe
  w.iwidth(2);                // Stellen in Indexausgabe
  w.width(20);                // Anzahl Stellen gesamt
  w.precision(16);            // Anzahl Nachkommastellen
  std::cout << w << std::endl;// nun mit mehr Stellen

  // Hilfsfunktionen
  zero(w);                    // das selbe wie w=0.0
  fill(w,(float)1.0);         // das selbe wie w=1.0
  fill(w,(float)0.0,(float)0.1); // w[0]=0, w[1]=0.1, w[2]=0.2, ...
  unitvector(w,2);            // kartesischer Einheitsvektor
  gnuplot("test.dat",w);      // gnuplot Ausgabe: i w[i]
  gnuplot("test2.dat",w,z);   // gnuplot Ausgabe: w[i] z[i]

  // Funktionsaufruf
  product(x);
  std::cout << "x=" << x << std::endl;
  std::cout << sum(x) << std::endl;
}
