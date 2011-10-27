#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

#include "hodgkinhuxley.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef HodgkinHuxley<Number> Model; // Model type
  Model model;                         // instantiate model

  typedef RKF45<Model> Solver;         // Solver type
  Solver solver(model);                // instantiate solver
  solver.set_dt(1.0/16.0);             // set initial time step
  solver.set_TOL(0.001);

  std::vector<Number> times;           // store time values here
  std::vector<Vector<Number> > states; // store states here
  std::vector<Number> dts;             // store delta t
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state
  dts.push_back(solver.get_dt());      // initial dt

  Number T = 130.0;
  while (solver.get_time()<T-1e-6) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      dts.push_back(solver.get_dt());      // used dt
    }

  gnuplot("hodgkinhuxley.dat",times,states,dts); // output model result

  return 0;
}
