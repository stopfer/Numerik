#include <iostream>
#include <vector>
#include "hdnum.hh"

using namespace hdnum;

#include "twobody.hh"
#include "restricted3BP.hh"
#include "threebody.hh"
#include "figureeight.hh"
#include "expliciteuler.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef TwoBody<Number> Model; // Model type
  //typedef ThreeBody<Number> Model; // Model type
  //typedef Restricted3Body<Number> Model; // Model type
  //typedef FigureEight<Number> Model; // Model type
  Model model;                         // instantiate model

  typedef RungeKutta4<Model> SubSolver;         // Solver type
  SubSolver subsolver(model);                // instantiate solver
  typedef RE<Model,SubSolver> Solver;         // Solver type
  Solver solver(model,subsolver);                // instantiate solver
  solver.set_dt(1.0/512.0);             // set initial time step
  solver.set_TOL(1E-10);

  std::vector<Number> times;           // store time values here
  std::vector<Vector<Number> > states; // store states here
  std::vector<Number> dts;             // store delta t
  times.push_back(solver.get_time());  // initial time
  states.push_back(solver.get_state()); // initial state
  dts.push_back(solver.get_dt());      // initial dt
  Number e_0(model.energy(solver.get_state()));
  std::cout << "initial energy: "  << std::scientific << std::showpoint 
	    << std::setprecision(12) << e_0 << std::endl;

  Number T = 100.0; // 100.0, figureeight: 2.1
  while (solver.get_time()<T-1e-8) // the time loop
    {
      solver.step();                  // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      dts.push_back(solver.get_dt());      // used dt
    }

  Number e_N(model.energy(solver.get_state()));
  std::cout << "final energy  : " << std::scientific << std::showpoint 
	    << std::setprecision(12) << e_N
	    << " |e_0-e_N|/|e_0|: " << std::scientific << std::showpoint 
	    << std::setprecision(12) << fabs(e_0-e_N)/fabs(e_0) << std::endl;
  std::cout << "number of f evaluations: " << model.get_count() << std::endl;
  gnuplot("nbody.dat",times,states,dts); // output model result
  //gnuplot("twobody.dat",times,states); // output model result
  //gnuplot("restricted3BP.dat",times,states); // output model result
  //gnuplot("threebody.dat",times,states); // output model result
  //gnuplot("figureeight.dat",times,states); // output model result

  return 0;
}
