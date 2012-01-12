#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include <gmpxx.h>
#include "hdnum.hh"

using namespace hdnum;

#include "model.hh"
#include "pc_incomplete.hh"

int main ()
{
  typedef double Number;               // define a number type

  typedef ConvergenceModel<Number> Model;  // Model type
  Model model(200.0,1./901. ,-3.);             // instantiate model

  Number h0 = 1/(pow(2,2));
  
  Vector<Number> times;           // store time values here
  Vector<Vector<Number> > states; // store states here

  for(int i=0; i<=5; i++){
    Number h = h0/pow(2,i);
    typedef RungeKutta4<Model> RKSolver; // RungeKutta solver
    RKSolver rksolver(model);                // instantiate solver
    rksolver.set_dt(h);
    times.push_back(rksolver.get_time());  // initial time
    states.push_back(rksolver.get_state()); // initial state

    typedef PC<Model> PCSolver; // Prädiktor/Korrektor solver

    Number t = -3;
    //do 3 RK steps to get inital values for LMM
    rksolver.step();
    times.push_back(rksolver.get_time());
    states.push_back(rksolver.get_state());
    Vector<Number> un3 = rksolver.get_state();
    rksolver.step();
    times.push_back(rksolver.get_time());
    states.push_back(rksolver.get_state());
    Vector<Number> un2 = rksolver.get_state();
    rksolver.step();
    times.push_back(rksolver.get_time());
    states.push_back(rksolver.get_state());
    Vector<Number> un1 = rksolver.get_state();
    PCSolver pcsolver(model,un3,un2,un1,h);    // instantiate PC solver with 
                                            //inital value from RK

    t = -3+3*h;
    while(t < 3.0){
      pcsolver.step(); //make time step
      times.push_back(pcsolver.get_time());
      states.push_back(pcsolver.get_state());
      t = pcsolver.get_time();
    }
    Number ut = 1.0 / (1.0+100.0*pcsolver.get_time()*pcsolver.get_time());
    std::cout << "h " << h << " eh_ab " << fabs(double(ut-(pcsolver.get_state()[0]))) << std::endl;

    char name[255];
    sprintf(name,"%s-%d","uebung9",i);
    gnuplot(name,times,states); // output model result
    times.clear();
    states.clear();
  }

  return 0;
}
