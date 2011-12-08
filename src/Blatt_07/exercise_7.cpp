#include <iostream>
#include <math.h>

#include "hdnum.hh"
#include "trapezoid_rule.h"

#include "initial_value_problem.h"


int main() {
   typedef double Number; // define a number type

   const Number t0 = 0.0; // initial time
   Number tStep = 1.0; // delta t
   const Number tMax = 2.0; // end time
   hdnum::Vector<Number> u0(3); // initial state
   hdnum::Vector<Number> uError(3); // initial state
   hdnum::Vector<Number> uExpected(3); // initial state
   hdnum::DenseMatrix<Number> A(3, 3, 0.0);
   const Number T = tMax-t0;
   const Number maxError = pow(10.0, -10.0);

   //set values of u0 and A
   u0[0] = 1.0; u0[1] = 0.0; u0[2] = -1.0;
   uExpected[0] = ((1.0/2.0) * exp(-2.0 * tMax)) + ((1.0/2.0) * exp(-40.0 * tMax) * (cos(40 * tMax) + sin(40 * tMax)));
   uExpected[1] = ((1.0/2.0) * exp(-2.0 * tMax)) - ((1.0/2.0) * exp(-40.0 * tMax) * (cos(40 * tMax) + sin(40 * tMax)));
   uExpected[2] = -exp(-40.0 * tMax) * (cos(40 * tMax) - sin(40 * tMax));
   A[0][0] = -21.0; A[0][1] = 19.0;  A[0][2] = -20.0;
   A[1][0] = 19.0;  A[1][1] = -21.0; A[1][2] = 20.0;
   A[2][0] = 40.0;  A[2][1] = -40.0; A[2][2] = -40.0;

   typedef InitialValueProblem<hdnum::DenseMatrix,hdnum::Vector, Number> Model; // Model type
   Model model(A, u0, t0); // instantiate model
   typedef Trapezoid<Model> Solver; // solver
   Solver solver(model); // instantiate solver
   int i = 0;
   do {
     tStep = pow(2.0, -i);
     std::cout << "Start new iteration: " << i << " tStep: " << tStep << std::endl;
     // initialize solver
     solver.set_dt(tStep);
     solver.set_state(t0, u0);

     while (solver.get_time() <= tMax) // the time loop
     {
       solver.step(); // advance model by one time step
     }
     uError = solver.get_state() - uExpected;
     i++;
   } while(uError.two_norm() > maxError);

   std::cout << "last state: " << solver.get_state() << std::endl;
   std::cout << "expected state: " << uExpected << std::endl;
   std::cout << "computed delta t: " << tStep * 2.0 << " = 2^(-" << i - 1 << ")" << std::endl;
   std::cout << "computed error (Euclidean norm): " << uError.two_norm() << std::endl;

   return 0;
}
