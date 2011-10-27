// funktion.cc
#include <iostream>

double f (double x)
{
  return x*x;
}

int main ()
{
  double x(2.0);
  std::cout << "f(" << x << ")=" << f(x) << std::endl;
}
