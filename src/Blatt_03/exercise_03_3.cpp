#include <iostream>

#include "hdnum.hh"
#include "initial_value_problem.h"
#include "expliciteuler.hh"             //from examplecode of hdnum
#include "heuns_method.h"

// y - eps < x < y + eps ?
bool inRange(const double x, const double y, const double eps) {
  if(x > y - eps) {
    if (x < y + eps) {
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}

template <class T>
T computeConvergenceOrder(const T x1, const T x2, const T xExact) {
  const T error1 = x1 - xExact;
  const T error2 = x2 - xExact;
  return log(fabs(error1 / error2)) / log(2.0);
}

template <class VECTOR>
typename VECTOR::value_type computeAverage(const VECTOR& vec) {
  typename VECTOR::value_type average = 0.0;
  for(size_t i = 0; i < vec.size(); i++) {
    average += vec.at(i);
  }
  average /= vec.size();
  return average;
}

int main() {
  typedef double Number;                    // define a number type

  const Number t0             = -0.5;       // initial time
        Number tStep          = 0.00001;    // delta t
  const Number tIntermediate  = 0.0;        // state at this time will be printed on console
  const Number uIntermediate  = sqrt(1 - (tIntermediate * tIntermediate)); // exact state at tIntermediate
  const Number tMax           = 1.0;        // end time
  const Number u0             = 0.75;       // initial state
  const Number lambda         = -1.0;     // factor for the model

  const size_t iStart         = 3;
  const size_t iEnd           = 8;

  hdnum::Vector<Number> statesIntermediate;     // store intermediate state values
  for(size_t i = iStart; i <= iEnd; i++) {
    tStep = pow(2.0, -static_cast<double>(i));
    typedef InitialValueProblem<Number> Model;  // Model type
    Model model(lambda, u0, t0);                // instantiate model

    typedef ExplicitEuler<Model> Solver;        // Solver type
    //typedef HeunsMethod<Model> Solver;        // Solver type
    Solver solver(model);                       // instantiate solver
    solver.set_dt(tStep);                       // set initial time step

    hdnum::Vector<Number> times;                // store time values here
    hdnum::Vector<hdnum::Vector<Number> > states;   // store states here
    times.push_back(solver.get_time());       // initial time
    states.push_back(solver.get_state());       // initial state

    while (solver.get_time() < tMax)         // the time loop
    {
      solver.step();                          // advance model by one time step
      times.push_back(solver.get_time());       // save time
      states.push_back(solver.get_state());     // and state
      if(inRange(tIntermediate, solver.get_time(), tStep/2.0)) {    // print state at intermediate time
        //std::cout << "State at time t(delta t = " << tStep << ") = " << solver.get_time() << ": " << solver.get_state() << std::endl;
        statesIntermediate.push_back(solver.get_state().at(0));
      }
    }
    if(i == iEnd) {
      gnuplot("exercise_03_3.dat",times,states);    // output model result
    }
  }

  hdnum::Vector<Number> convergenceOrders;     // store computed convergence orders
  for(size_t i = 0; i < iEnd - iStart; i++) {
    convergenceOrders.push_back(computeConvergenceOrder(statesIntermediate.at(i), statesIntermediate.at(i + 1), uIntermediate));
    std::cout << "convergence order for h = " << pow(2.0, -static_cast<double>(i + iStart)) << " : " << convergenceOrders.back() << std::endl;
  }

  Number averageConvergenceOrder = computeAverage(convergenceOrders);
  std::cout << "average convergence order: " << averageConvergenceOrder << std::endl;

  return 0;
}
