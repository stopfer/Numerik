#include <iostream>
#include <vector>
#include "hdnum.hh"

#include "modelproblem.hh"
#include "expliciteuler.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef ModelProblem<Number> Model;  // Model type
  Model model(-100.0);                   // instantiate model

  typedef ExplicitEuler<Model> Solver; // Solver type
  Solver solver(model);                // instantiate solver
  solver.set_dt(0.02);                  // set initial time step

  hdnum::Vector<Number> times;           // store time values here
  hdnum::Vector<hdnum::Vector<Number> > states; // store states here
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state

  while (solver.get_time()<5.0-1e-6) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
    }

  gnuplot("mp2-ee-0.02.dat",times,states); // output model result

  return 0;
}
