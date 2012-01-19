#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include "gmp.h"
#include "hdnum.hh"

using namespace hdnum;

#include "model.hh"
#include "adamsbashfort.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef ConvergenceModel<Number> Model;  // Model type
  Model model(200.0,1./901. ,-3.);             // instantiate model

  Number h0 = 1/(pow(2,2));
  
  for(int i=0; i<=6; i++){
    model.counter = 0;
    Number h = h0/pow(2,i);
    typedef RungeKutta4<Model> RKSolver; // RungeKutta solver
    RKSolver rksolver(model);                // instantiate rksolver for absolver (first steps)
    rksolver.set_dt(h);
	RKSolver rksolver2(model);                // instantiate rksolver2 as comparison
	rksolver2.set_dt(h);  
    typedef AdamsBashfort<Model> ABSolver; // AdamsBashfort solver

    Number t = -3;
    rksolver.step();
    Vector<Number> u2 = rksolver.get_state();
    rksolver.step();
    Vector<Number> u3 = rksolver.get_state();
    rksolver.step();
    Vector<Number> u4 = rksolver.get_state();
    ABSolver absolver(model,u2,u3,u4,h);                // instantiate absolver

	Number abfvalues = 12;
    t = -3+3*h;
    while(t<0){
      absolver.step(); //make time step with AdamsBashforth
	  abfvalues++;
      t = absolver.get_time();
    }
	
	Number rkfvalues = 0;
	Number t2 = -3;
	while(t2<0){
	  rksolver2.step(); //make time step with RungeKutta4
	  rkfvalues += 4;
	  t2 = rksolver2.get_time();
	}
	  
    Number ut = 1;
    std::cout << "h " << h << " eh_ab " << fabs(ut-(absolver.get_state())[0]) << std::endl;
	std::cout << "abfvalues " << abfvalues << std::endl;
    std::cout << "h " << h << " eh_rk " << fabs(ut-(rksolver2.get_state())[0]) << std::endl; 	
	std::cout << "rkfvalues " << rkfvalues << std::endl;
  }

  return 0;
}
