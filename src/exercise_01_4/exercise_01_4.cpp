#include <iostream>

#include "hdnum.hh"
#include "initial_value_problem.h"
#include "expliciteuler.hh"

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

int main() {
	typedef double Number;               			// define a number type

	const Number t0 = -3.0;							// initial time
	const Number tStep = 0.00001;						// delta t
	const Number tIntermediate = 1.0;				// state at this time will be printed on console
	const Number tMax = 3.0;						// end time
	const Number u0 = 1.0/901.0;					// initial state
	const Number lambda = -200.0;					// factor for the model

	typedef InitialValueProblem<Number> Model;  	// Model type
	Model model(lambda, u0, t0);					// instantiate model

	typedef ExplicitEuler<Model> Solver; 			// Solver type
	Solver solver(model);                			// instantiate solver
	solver.set_dt(tStep);                  			// set initial time step

	hdnum::Vector<Number> times;           			// store time values here
	hdnum::Vector<hdnum::Vector<Number> > states; 	// store states here
	times.push_back(solver.get_time());  			// initial time
	states.push_back(solver.get_state()); 			// initial state

	while (solver.get_time() <= tMax) 				// the time loop
	{
	  solver.step();                  				// advance model by one time step
	  times.push_back(solver.get_time()); 			// save time
	  states.push_back(solver.get_state()); 		// and state
	  if(inRange(tIntermediate, solver.get_time(), tStep/2.0)) {		// print state at intermediate time
		  std::cout << "State at time t = " << solver.get_time() << ": " << solver.get_state() << std::endl;
	  }
	}

	gnuplot("exercise_01_4.dat",times,states); 		// output model result

	return 0;
}
