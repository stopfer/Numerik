#include <iostream>
#include <vector>
#include <math.h>
#include <gmpxx.h>
#include "hdnum.hh"

using namespace hdnum;

#include "modelproblem.hh"
#include "expliciteuler.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(3.0);                   // instantiate model

  typedef Heun3<Model> Solver;   // Solver type

  Number dt=0.25;
  Number last=0.0;

  for (int i=0; i<15; i++)
    { 
      Solver solver(model);                // instantiate solver
      solver.set_dt(dt);                  // set initial time step

      while (solver.get_time()<1.0-1e-8) // the time loop
	solver.step();                   // advance model by one time step

      Number error = fabs(exp(3.0)-solver.get_state()[0]);

      std::cout << std::setw(2) << i 
		<< " dt=" << std::scientific << std::showpoint << std::setprecision(10) << dt 
		<< " error=" << std::scientific << std::showpoint 
		<< std::setprecision(10) << error 
		<< " reduction=" << std::scientific << std::showpoint 
		<< std::setprecision(10) << last/error
		<< std::endl;
      last = error;
      dt *= 0.5;
    }

  return 0;
}
