#include <iostream>
#include <gmpxx.h>
#include "hdnum.hh"

using namespace hdnum;

const int n=7;
typedef mpf_class number;

int main ()
{
  mpf_set_default_prec(1024);

  Vector<number> x(n);
  Vector<number> b(n);
  Vector<number> s(n);
  Array<std::size_t> p(n);
  Array<std::size_t> q(n);
  DenseMatrix<number> A(n,n);
  fill(x,number(1.0),number(1.0));
  vandermonde(A,x);
  A.mv(b,x);
  row_equilibrate(A,s);
  x = number(0.0);
  lr_fullpivot(A,p,q);
  apply_equilibrate(s,b);
  permute_forward(p,b);
  solveL(A,b,b);
  solveR(A,x,b);
  permute_backward(q,x);
  std::cout << x << std::endl;
}
