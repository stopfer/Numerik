// matrizen.cc 
#include <iostream>    // notwendig zur Ausgabe
#include "hdnum.hh"    // hdnum header

// Beispiel wie man A und b für ein
// Gleichungssystem initialisieren könnte
template<class T>
void initialize (hdnum::DenseMatrix<T> &A, hdnum::Vector<T> &b)
{
  if (A.rowsize()!=A.colsize() || A.rowsize()==0) 
    HDNUM_ERROR("need square and nonempty matrix");
  if (A.rowsize()!=b.size())
    HDNUM_ERROR("b must have same size as A");
  for (int i=0; i<A.rowsize(); ++i)
  {
    b[i] = 1.0;
      for (int j=0; j<A.colsize(); ++j)
    if (j<=i) A[i][j]=1.0; else A[i][j]=0.0;
  }
}

int main ()
{
  // Konstruktion
  hdnum::DenseMatrix<float> A;            // leere Matrix mit Größe 0x0
  hdnum::DenseMatrix<float> B(10,10);     // 10x10 Matrix uninitialisiert
  hdnum::DenseMatrix<float> C(10,10,0.0); // 10x10 Matrix initialisiert

  // Zugriff auf Vektorelemente
  for (int i=0; i<B.rowsize(); ++i)
    for (int j=0; j<B.colsize(); ++j)
      B[i][j] = 0.0;          // jetzt ist B initialisiert

  // Kopie und Zuweisung
  hdnum::DenseMatrix<float> D(B);    // D ist eine Kopie von B
  A = D;                      // A ist nun identisch mit D!
  A[0][0] = 3.14;
  B[0][0] = 3.14;
  
  // Rechnen mit Matrizen und Vektoren
  A += B;                     // A = A+B
  A -= B;                     // A = A-B
  A *= 1.23;                  // Multiplikation mit Skalar
  A /= 1.23;                  // Division durch Skalar
  A.update(1.23,B);           // A = A + s*B

  hdnum::Vector<float> x(10,1.0);    // make two vectors
  hdnum::Vector<float> y(10,2.0);
  A.mv(y,x);                  // y = A*x
  A.umv(y,x);                 // y = y + A*x
  A.umv(y,(float)-1.0,x);     // y = y + s*A*x
  C.mm(A,B);                  // C = A*B
  C.umm(A,B);                 // C = C + A*B

  // Ausgabe
  A.iwidth(2);                // Stellen in Indexausgabe
  A.width(11);                // Anzahl Stellen gesamt
  A.precision(4);             // Anzahl Nachkommastellen
  std::cout << A << std::endl;// schöne Ausgabe

  // Hilfsfunktionen
  identity(A);                // setze A auf Einheitsmatrix
  std::cout << A << std::endl;
  spd(A);                     // eine s.p.d. Matrix
  std::cout << A << std::endl;
  fill(x,(float)1,(float)1);
  vandermonde(A,x);           // Vandermondematrix
  std::cout << A << std::endl;
}
