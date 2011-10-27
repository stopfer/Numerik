#include <iostream>
#include <iomanip>

int main()
{
  int n;
  std::cout << "Geben Sie eine natürliche Zahl (kleiner als 46) ein: ";
  std::cin >> n;

  int fib = 0;

  if( n == 0 )
  {
    fib = 0;
    std::cout << "fib(0) = " << fib << std::endl;
  }
  else if( n == 1 )
  {
    fib = 1;
    std::cout << "fib(0) = " << 0 << std::endl;
    std::cout << "fib(1) = " << fib << std::endl;
  }
  else
  {
    int f0 = 0;     // Startwert 1
    int f1 = 1;     // Startwert 2
    fib = f1 + f0;
    int k=2;
    
    std::cout << "fib( 0) = " << std::setw(10) << f0 << std::endl;
    std::cout << "fib( 1) = " << std::setw(10) << f1 << std::endl;
    while(k<=n)
    {
      f0 = f1;           // Vor-Vorgaenger
      f1 = fib;          // Vorgaenger
      fib = f1 + f0;
      std::cout << "fib(" << std::setw(2)  <<  k << ") = " 
          << std::setw(10) << fib << std::endl;
      k=k+1;
    }
  }
}
