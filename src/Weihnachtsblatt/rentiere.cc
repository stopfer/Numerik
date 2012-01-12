#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include <gmpxx.h>
#include "hdnum.hh"

#include "rentiere_model.hh"

using namespace hdnum;

int main () {
  // define a number type
  typedef double Number;
  // define constants
  const Number t0 = 0.0;
  const Number tMax = 358.0;
  const Number dt = 1.0;
  const Number u0 = 100.0;

  // model type
  typedef ReindeerModel<Number> Model;
  // instantiate model
  Model model(u0, t0);

  // storage for computed times and states
  Vector<Number> times;
  Vector<Vector<Number> > states;

  // solver type
  typedef RungeKutta4<Model> RKSolver;
  // instantiate solver
  RKSolver rksolver(model);
  rksolver.set_dt(dt);

  // set initial time and state
  times.push_back(rksolver.get_time());
  states.push_back(rksolver.get_state());

  // compute times and states
  Number t = t0;
  while(t < tMax){
    rksolver.step(); //make time step
    times.push_back(rksolver.get_time());
    states.push_back(rksolver.get_state());
    t = rksolver.get_time();
  }

  // write computed results in file
  gnuplot("reindeerpopulation.dat",times,states); // output model result
  std::cout << "reindeer population at time t = " << times.back() << ": " << states.back() << std::endl;
  return 0;
}
