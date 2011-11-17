#include <iostream>

#include "hdnum.hh"
#include "initial_value_problem.h"

int main() {
   typedef double Number; // define a number type

   const Number t0 = -3.0; // initial time
   const Number tStep0 = 0.0001; // delta t
   const Number tMax = -1.0; // end time
   const Number u0 = 1.0 / 901.0; // initial state
   const Number T= tMax-t0;

   typedef InitialValueProblem<Number> Model; // Model type
   Model model( u0, t0); // instantiate model
   typedef hdnum::Kutta3<Model> Solver; // solver
   Solver solver(model); // instantiate solver


   hdnum::Vector<Number> times; // store time values here
   hdnum::Vector<hdnum::Vector<Number> > states; // store states here
   times.push_back(solver.get_time()); // initial time
   states.push_back(solver.get_state()); // initial state
   Number h_n=tStep0;
   Number maxError=0;
   while (solver.get_time() <= tMax) // the time loop
   {
      std::cout<<solver.get_time()<<"\n";
      //do 2 steps with half stepsize
      solver.set_dt(h_n/2.0); // set initial time step
      solver.step(); // advance model by a half time step
      solver.step(); // advance model by a half time step
      Number yHalfStep=solver.get_state()[0];
      // 
      solver.set_dt(h_n); // set normal time step
      solver.set_state(times.back(),states.back());
      solver.step(); // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state

      Number yFullStep=solver.get_state()[0];
      Number absDiffHalfFullStep=std::fabs(yHalfStep-yFullStep);

      const size_t m=3; //kutta 3
      Number epsilon=std::pow(10.0,-10.0);
      Number TOL=epsilon*std::fabs(solver.get_state()[0])/h_n;
      Number TermA=std::pow(2.0*h_n,m+1)*(1.0-std::pow(2.0,-1.0*m));
      Number TermB=T*std::fabs(absDiffHalfFullStep);
      Number h_opt=((TermA*TOL)/TermB);
      h_opt=std::pow(h_opt, Number(1)/Number(m));
      
      //h_n=h_opt  is fine with  1/2 h_n <=h_opt
      //do the "real" step with the optimized h_n
      h_n=h_opt;
      solver.set_dt(h_n);
      solver.step(); // advance model by one time step
      times.push_back(solver.get_time()); // save time
      states.push_back(solver.get_state()); // and state
      const Number error=fabs(1.0/(1.0+100.0*std::pow(solver.get_time(),2))-solver.get_state()[0]);
      if(error>maxError){
         maxError=error;
      }
   }
   
   std::cout<<"max error: "<<maxError<<"\n";
   solver.set_state(t0,states.front());
   solver.set_dt(tStep0);
   h_n=tStep0;
   while (solver.get_time() <= tMax) // the time loop
   {
      std::cout<<solver.get_time()<<"\n";
      solver.set_dt(h_n);
      solver.step();
      const Number error=fabs(1.0/(1.0+100.0*std::pow(solver.get_time(),2))-solver.get_state()[0]);
      if(error>maxError){
        //std::cout<<" error h_old: "<<h_n<<"h_new"<<h_n/2<<"\n";
        h_n*=0.99;
        solver.set_dt(h_n);
        solver.set_state(t0,states.front());
      }
   }
   std::cout<<"fixed h_n for same error : "<<h_n<<"\n";

   return 0;
}
