#include <assert.h>
#include "hdnum.hh"   

using namespace hdnum; // Namen ohne hdnum:: verwenden 

template<class M>
class PC
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

public:

  PC (const M& model_, Vector<number_type> un3_, 
		 Vector<number_type> un2_, Vector<number_type> un1_,
		 number_type dt_)
    : model(model_), u(model.size()), un0(model.size()), un1(model.size()),
      un2(model.size()), un3(model.size()), un4(model.size()),
      fn0(model.size()), fn1(model.size()), fn2(model.size()),
      fn3(model.size()), fn4(model.size())
  {
    model.initialize(tn1,un4); // time is currently t=t0, determine u_(n-4)
    un3 = un3_;                // u_(n-3)
    un2 = un2_;                // u_(n-2)
    un1 = un1_;                // u_(n-1)
    dt = dt_;                  // timestep
    t0 = tn1;                    // save starting time
    tn1 = tn1+3*dt;              // we need 3 RK steps first,
                               // so we start with tn-1=t0+3*dt
  }

  //! do one step
  void step () {
    //Prädiktor
    if ( tn1==t0+3*dt ) {      // first step, so we need to compute all fs
      model.f(tn1-3*dt,un4,fn4);  // f_(n-4)(t_(n-4))
      model.f(tn1-2*dt,un3,fn3);  // f_(n-3)(t_(n-3))
      model.f(tn1-1*dt,un2,fn2);  // f_(n-2)(t_(n-2))
      model.f(tn1,un1,fn1);       // f_(n-1)(t_(n-1))
    }
    
    /************
     ADD MISSING CODE HERE:
     compute u with 4th order Adams-Bashfort:
     u = u_(n-1) + 1/24*h*(55f_(n-1)-59f_(n-2)+37f_(n-3)-9f_(n-4))
     
     useful method:
     u.update(a,v) => u = u + a*v (a scalar, u,v vector)
    *************/
    // compute intermediate value to update u
    Vector<number_type> fnsum(model.size(), 0.0);
    fnsum.update(55.0, fn1);
    fnsum.update(-59.0, fn2);
    fnsum.update(37.0, fn3);
    fnsum.update(-9.0, fn4);
    // update u
    u = un1;
    u.update(((1.0/24.0) * dt), fnsum);

    // Corrector
    Vector<number_type> fnstar(model.size());

    /************
     ADD MISSING CODE HERE: compute u with 4th order Adams-Moulton: 
     u= u_(n-1) + 1/24*h*(9fnstar+19f_(n-1)-5f_(n-2)+f_(n-3))
     with fnstar=f(t_n,u), where u is the old solution from Prädiktor
     Adams-Bashfort.  

     Then compute new fnstar=f(t_n,u) with the new u from Corrector

    *************/
    // compute fnstar
    model.f(tn1 + dt, u, fnstar);
    // reset fnsum
    fnsum = 0.0;
    // compute intermediate values to update u
    fnsum.update(9.0, fnstar);
    fnsum.update(19.0, fn1);
    fnsum.update(-5.0, fn2);
    fnsum.update(1.0, fn3);
    // update u
    u = un1;
    u.update(((1.0/24.0) * dt), fnsum);
    // compute new fnstar
    model.f(tn1 + dt, u, fnstar);

    // prepare for next step: (n-x) -> (n-(x+1)) 
    un4=un3;
    un3=un2;
    un2=un1;
    un1=u;
    fn4=fn3;
    fn3=fn2;
    fn2=fn1;
    fn1=fnstar;
    tn1 += dt;
    
  }

  //! get current state
  const Vector<number_type>& get_state () const
  {
    return u;
  }

  //! get current time
  time_type get_time () const
  {
    return tn1;
  }

  //! get dt used in last step (i.e. to compute current state)
  time_type get_dt () const
  {
    return dt;
  }

private:
  const M& model;                    
  time_type tn1, t0, dt;
  Vector<number_type> u;
  Vector<number_type> un0;
  Vector<number_type> un1;
  Vector<number_type> un2;
  Vector<number_type> un3;
  Vector<number_type> un4;
  Vector<number_type> fn0;
  Vector<number_type> fn1;
  Vector<number_type> fn2;
  Vector<number_type> fn3;
  Vector<number_type> fn4;
};
