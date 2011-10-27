// string.cc
#include <iostream>
#include <string>
int main ()
{
  std::string m1("Zeichen");
  std::string leer("   ");
  std::string m2("kette");
  std::cout << m1+leer+m2 << std::endl;
}
