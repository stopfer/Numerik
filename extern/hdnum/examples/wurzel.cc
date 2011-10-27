#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

template<class N>
class WurzelProblem
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  WurzelProblem (number_type a_)
    : a(a_)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 1;
  }

  //! model evaluation
  void F (const Vector<N>& x, Vector<N>& result) const
  {
    result[0] = x[0]*x[0] - a;
  }

  //! jacobian evaluation needed for implicit solvers
  void F_x (const Vector<N>& x, DenseMatrix<N>& result) const
  {
    result[0][0] = number_type(2.0)*x[0];
  }

private:
  number_type a;
};



int main ()
{
  typedef double Number;                 // Zahlentyp

  typedef WurzelProblem<Number> Problem; // Problemtyp
  Problem problem(2.0);                  // Eine Instanz des Problems

  Newton newton;                         // Ein Newtonobjekt
  newton.set_maxit(20);                  // Setze diverse Parameter
  newton.set_verbosity(2);    
  newton.set_reduction(1e-12);
  newton.set_abslimit(1e-15);
  newton.set_linesearchsteps(3);  

  Vector<Number> u(problem.size());      // Objekt für die Loesung
  u[0] = 17.0;                           // Startwert
  newton.solve(problem,u);               // Berechne Lösung
  std::cout << "Ergebnis: " << u[0] << std::endl;

  return 0;
}
