// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef HDNUM_ODE_HH
#define HDNUM_ODE_HH

#include<vector>
#include "newton.hh"

/** @file
 *  @brief solvers for ordinary differential equations
 */

namespace hdnum {

  /** @brief Explicit Euler method as an example for an ODE solver

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class EE
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    EE (const M& model_)
      : model(model_), u(model.size()), f(model.size())
    {
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      model.f(t,u,f);   // evaluate model
      u.update(dt,f);   // advance state
      t += dt;          // advance time
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 1;
    }
    
  private:
    const M& model;
    time_type t, dt;
    Vector<number_type> u;
    Vector<number_type> f;
  };

  /** @brief Modified Euler method (order 2 with 2 stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class ModifiedEuler
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    ModifiedEuler (const M& model_)
      : model(model_), u(model.size()), w(model.size()), k1(model.size()), k2(model.size())
    {
      c2 = 0.5;
      a21 = 0.5;
      b2 = 1.0;
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // final
      u.update(dt*b2,k2);
      t += dt;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 2;
    }
    
  private:
    const M& model;
    time_type t, dt;
    time_type c2,a21,b2;
    Vector<number_type> u,w;
    Vector<number_type> k1,k2;
  };


  /** @brief Heun method (order 2 with 2 stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class Heun2
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    Heun2 (const M& model_)
      : model(model_), u(model.size()), w(model.size()), k1(model.size()), k2(model.size())
    {
      c2 = 1.0;
      a21 = 1.0;
      b1 = 0.5;
      b2 = 0.5;
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // final
      u.update(dt*b1,k1);
      u.update(dt*b2,k2);
      t += dt;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 2;
    }
    
  private:
    const M& model;
    time_type t, dt;
    time_type c2,a21,b1,b2;
    Vector<number_type> u,w;
    Vector<number_type> k1,k2;
  };


  /** @brief Heun method (order 3 with 3 stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class Heun3
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    Heun3 (const M& model_)
      : model(model_), u(model.size()), w(model.size()), k1(model.size()), 
        k2(model.size()), k3(model.size())
    {
      c2 = time_type(1.0)/time_type(3.0);
      c3 = time_type(2.0)/time_type(3.0);
      a21 = time_type(1.0)/time_type(3.0);
      a32 = time_type(2.0)/time_type(3.0);
      b1 = 0.25;
      b2 = 0.0;
      b3 = 0.75;
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // stage 3
      w = u;
      w.update(dt*a32,k2);
      model.f(t+c3*dt,w,k3);

      // final
      u.update(dt*b1,k1);
      u.update(dt*b3,k3);
      t += dt;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 3;
    }
    
  private:
    const M& model;
    time_type t, dt;
    time_type c2,c3,a21,a31,a32,b1,b2,b3;
    Vector<number_type> u,w;
    Vector<number_type> k1,k2,k3;
  };

  /** @brief Kutta method (order 3 with 3 stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class Kutta3
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    Kutta3 (const M& model_)
      : model(model_), u(model.size()), w(model.size()), k1(model.size()), 
        k2(model.size()), k3(model.size())
    {
      c2 = 0.5;
      c3 = 1.0;
      a21 = 0.5;
      a31 = -1.0;
      a32 = 2.0;
      b1 = time_type(1.0)/time_type(6.0);
      b2 = time_type(4.0)/time_type(6.0);
      b3 = time_type(1.0)/time_type(6.0);
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // stage 3
      w = u;
      w.update(dt*a31,k1);
      w.update(dt*a32,k2);
      model.f(t+c3*dt,w,k3);

      // final
      u.update(dt*b1,k1);
      u.update(dt*b2,k2);
      u.update(dt*b3,k3);
      t += dt;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 3;
    }
    
  private:
    const M& model;
    time_type t, dt;
    time_type c2,c3,a21,a31,a32,b1,b2,b3;
    Vector<number_type> u,w;
    Vector<number_type> k1,k2,k3;
  };

  /** @brief classical Runge-Kutta method (order 4 with 4 stages)

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class RungeKutta4
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    RungeKutta4 (const M& model_)
      : model(model_), u(model.size()), w(model.size()), k1(model.size()), 
        k2(model.size()), k3(model.size()), k4(model.size())
    {
      c2 = 0.5;
      c3 = 0.5;
      c4 = 1.0;
      a21 = 0.5;
      a32 = 0.5;
      a43 = 1.0;
      b1 = time_type(1.0)/time_type(6.0);
      b2 = time_type(2.0)/time_type(6.0);
      b3 = time_type(2.0)/time_type(6.0);
      b4 = time_type(1.0)/time_type(6.0);
      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // stage 3
      w = u;
      w.update(dt*a32,k2);
      model.f(t+c3*dt,w,k3);

      // stage 4
      w = u;
      w.update(dt*a43,k3);
      model.f(t+c4*dt,w,k4);

      // final
      u.update(dt*b1,k1);
      u.update(dt*b2,k2);
      u.update(dt*b3,k3);
      u.update(dt*b4,k4);
      t += dt;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 4;
    }
    
  private:
    const M& model;
    time_type t, dt;
    time_type c2,c3,c4,a21,a32,a43,b1,b2,b3,b4;
    Vector<number_type> u,w;
    Vector<number_type> k1,k2,k3,k4;
  };

  /** @brief Adaptive Runge-Kutta-Fehlberg method

      \tparam M the model type
  */
  template<class M>
  class RKF45
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    RKF45 (const M& model_)
      : model(model_), u(model.size()), w(model.size()), ww(model.size()), k1(model.size()), 
        k2(model.size()), k3(model.size()), k4(model.size()), k5(model.size()), k6(model.size()),
        steps(0), rejected(0)
    {
      TOL = time_type(0.0001);
      rho = time_type(0.8);
      alpha = time_type(0.25);
      beta = time_type(4.0);
      dt_min = 1E-12;

      c2 = time_type(1.0)/time_type(4.0);
      c3 = time_type(3.0)/time_type(8.0);
      c4 = time_type(12.0)/time_type(13.0);
      c5 = time_type(1.0);
      c6 = time_type(1.0)/time_type(2.0);

      a21 = time_type(1.0)/time_type(4.0);

      a31 = time_type(3.0)/time_type(32.0);
      a32 = time_type(9.0)/time_type(32.0);

      a41 = time_type(1932.0)/time_type(2197.0);
      a42 = time_type(-7200.0)/time_type(2197.0);
      a43 = time_type(7296.0)/time_type(2197.0);

      a51 = time_type(439.0)/time_type(216.0);
      a52 = time_type(-8.0);
      a53 = time_type(3680.0)/time_type(513.0);
      a54 = time_type(-845.0)/time_type(4104.0);

      a61 = time_type(-8.0)/time_type(27.0);
      a62 = time_type(2.0);
      a63 = time_type(-3544.0)/time_type(2565.0);
      a64 = time_type(1859.0)/time_type(4104.0);
      a65 = time_type(-11.0)/time_type(40.0);

      b1 = time_type(25.0)/time_type(216.0);
      b2 = time_type(0.0);
      b3 = time_type(1408.0)/time_type(2565.0);
      b4 = time_type(2197.0)/time_type(4104.0);
      b5 = time_type(-1.0)/time_type(5.0);

      bb1 = time_type(16.0)/time_type(135.0);
      bb2 = time_type(0.0);
      bb3 = time_type(6656.0)/time_type(12825.0);
      bb4 = time_type(28561.0)/time_type(56430.0);
      bb5 = time_type(-9.0)/time_type(50.0);
      bb6 = time_type(2.0)/time_type(55.0);

      model.initialize(t,u);
      dt = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! set tolerance for adaptive computation
    void set_TOL (time_type TOL_)
    {
      TOL = TOL_;
    }

    //! do one step
    void step ()
    {
      steps++;

      // stage 1
      model.f(t,u,k1);

      // stage 2
      w = u;
      w.update(dt*a21,k1);
      model.f(t+c2*dt,w,k2);

      // stage 3
      w = u;
      w.update(dt*a31,k1);
      w.update(dt*a32,k2);
      model.f(t+c3*dt,w,k3);

      // stage 4
      w = u;
      w.update(dt*a41,k1);
      w.update(dt*a42,k2);
      w.update(dt*a43,k3);
      model.f(t+c4*dt,w,k4);

      // stage 5
      w = u;
      w.update(dt*a51,k1);
      w.update(dt*a52,k2);
      w.update(dt*a53,k3);
      w.update(dt*a54,k4);
      model.f(t+c5*dt,w,k5);

      // stage 6
      w = u;
      w.update(dt*a61,k1);
      w.update(dt*a62,k2);
      w.update(dt*a63,k3);
      w.update(dt*a64,k4);
      w.update(dt*a65,k5);
      model.f(t+c6*dt,w,k6);

      // compute order 4 approximation
      w = u;
      w.update(dt*b1,k1);
      w.update(dt*b2,k2);
      w.update(dt*b3,k3);
      w.update(dt*b4,k4);
      w.update(dt*b5,k5);

      // compute order 5 approximation
      ww = u;
      ww.update(dt*bb1,k1);
      ww.update(dt*bb2,k2);
      ww.update(dt*bb3,k3);
      ww.update(dt*bb4,k4);
      ww.update(dt*bb5,k5);
      ww.update(dt*bb6,k6);

      // estimate local error
      w -= ww;
      time_type error(norm(w));
      time_type dt_opt(dt*pow(rho*TOL/error,0.2));
      dt_opt = std::min(beta*dt,std::max(alpha*dt,dt_opt));
      //std::cout << "est. error=" << error << " dt_opt=" << dt_opt << std::endl;

      if (error<=TOL)
        {
          t += dt;
          u = ww;
          dt = dt_opt;
        }
      else
        {
          rejected++;
          dt = dt_opt;
          if (dt>dt_min) step();
        }
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 5;
    }
    
    //! print some information
    void get_info () const
    {
      std::cout << "RE: steps=" << steps << " rejected=" << rejected << std::endl;
    }

  private:
    const M& model;
    time_type t, dt;
    time_type TOL,rho,alpha,beta,dt_min;
    time_type c2,c3,c4,c5,c6;
    time_type a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65;
    time_type b1,b2,b3,b4,b5; // 4th order
    time_type bb1,bb2,bb3,bb4,bb5,bb6; // 5th order
    Vector<number_type> u,w,ww;
    Vector<number_type> k1,k2,k3,k4,k5,k6;
    mutable size_type steps, rejected;
  };


  /** @brief Adaptive one-step method using Richardson extrapolation

      \tparam M a model
      \tparam S any of the (non-adaptive) one step methods (solving model M)
  */
  template<class M, class S>
  class RE
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    RE (const M& model_, S& solver_)
      : model(model_), solver(solver_), u(model.size()), 
        wlow(model.size()), whigh(model.size()), ww(model.size()),
        steps(0), rejected(0)
    {
      model.initialize(t,u); // initialize state
      dt = 0.1;              // set initial time step
      two_power_m = 1.0;
      for (size_type i=0; i<solver.get_order(); i++) 
        two_power_m *= 2.0;
      TOL = time_type(0.0001);
      rho = time_type(0.8);
      alpha = time_type(0.25);
      beta = time_type(4.0);
      dt_min = 1E-12;
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! set tolerance for adaptive computation
    void set_TOL (time_type TOL_)
    {
      TOL = TOL_;
    }

    //! do one step
    void step ()
    {
      // count steps done
      steps++;

      // do 1 step with 2*dt
      time_type H(2.0*dt);
      solver.set_state(t,u);
      solver.set_dt(H);
      solver.step();
      wlow = solver.get_state();

      // do 2 steps with dt
      solver.set_state(t,u);
      solver.set_dt(dt);
      solver.step();
      solver.step();
      whigh = solver.get_state();

      // estimate local error
      ww = wlow;
      ww -= whigh;
      time_type error(norm(ww)/(pow(H,1.0+solver.get_order())*(1.0-1.0/two_power_m)));
      time_type dt_opt(pow(rho*TOL/error,1.0/((time_type)solver.get_order())));
      dt_opt = std::min(beta*dt,std::max(alpha*dt,dt_opt));
      //std::cout << "est. error=" << error << " dt_opt=" << dt_opt << std::endl;

      if (dt<=dt_opt)
        {
          t += H;
          u = whigh;
          u *= two_power_m;
          u -= wlow;
          u /= two_power_m-1.0;
          dt = dt_opt;
        }
      else
        {
          rejected++;
          dt = dt_opt;
          if (dt>dt_min) step();
        }
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return solver.get_order()+1;
    }

    //! print some information
    void get_info () const
    {
      std::cout << "RE: steps=" << steps << " rejected=" << rejected << std::endl;
    }
    
  private:
    const M& model;
    S& solver;
    time_type t, dt;
    time_type two_power_m;
    Vector<number_type> u,wlow,whigh,ww;
    time_type TOL,rho,alpha,beta,dt_min;
    mutable size_type steps, rejected;
  };


  /** @brief Implicit Euler using Newton's method to solve nonlinear system

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
      \tparam S nonlinear solver
  */
  template<class M, class S>
  class IE
  {
    //! class providing nonlinear problem to be solved 
    // h_n f(t_n, y_n) - y_n + y_{n-1} = 0
    class NonlinearProblem
    {
    public:
      /** \brief export size_type */
      typedef typename M::size_type size_type;
      
      /** \brief export number_type */
      typedef typename M::number_type number_type;

      //! constructor stores parameter lambda
      NonlinearProblem (const M& model_, const Vector<number_type>& yold_, 
                        typename M::time_type tnew_, typename M::time_type dt_)
        : model(model_), yold(yold_), tnew(tnew_), dt(dt_)
      {}

      //! return number of componentes for the model
      std::size_t size () const
      {
        return model.size();
      }

      //! model evaluation
      void F (const Vector<number_type>& x, Vector<number_type>& result) const
      {
        model.f(tnew,x,result);
        result *= dt;
        result -= x;
        result += yold;
      }

      //! jacobian evaluation needed for implicit solvers
      void F_x (const Vector<number_type>& x, DenseMatrix<number_type>& result) const
      {
        model.f_x(tnew,x,result);
        result *= dt;
        for (size_type i=0; i<model.size(); i++) result[i][i] -= number_type(1.0);
      }

      void set_tnew_dt (typename M::time_type tnew_, typename M::time_type dt_)
      {
        tnew = tnew_;
        dt = dt_;
      }

    private:
      const M& model;
      const Vector<number_type>& yold;
      typename M::time_type tnew;
      typename M::time_type dt;
    };

  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    IE (const M& model_, const S& newton_)
      : verbosity(0), model(model_), newton(newton_), u(model.size()), unew(model.size())
    {
      model.initialize(t,u);
      dt = dtmax = 0.1; 
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dtmax = dt_;
    }

    //! set verbosity level
    void set_verbosity (size_type verbosity_)
    {
      verbosity = verbosity_;
    }

    //! do one step
    void step ()
    {
      if (verbosity>=2) 
        std::cout << "IE: step" << " t=" << t << " dt=" << dt << std::endl;
      NonlinearProblem nlp(model,u,t+dt,dt);
      bool reduced = false;
      error = false;
      while (1)
        {
          unew = u;
          newton.solve(nlp,unew);
          if (newton.has_converged())
            {
              u = unew;
              t += dt;
              if (!reduced && dt<dtmax-1e-13)
                {
                  dt = std::min(2.0*dt,dtmax);
                  if (verbosity>0) 
                    std::cout << "IE: increasing time step to " << dt << std::endl;
                }
              return;
            }
          else
            {
              if (dt<1e-12)
                {
                  HDNUM_ERROR("time step too small in implicit Euler");
                  error = true;
                  break;
                }
              dt *= 0.5;
              reduced = true;
              nlp.set_tnew_dt(t+dt,dt);
              if (verbosity>0) std::cout << "IE: reducing time step to " << dt << std::endl;
            }
        }
    }

    //! get current state
    bool get_error () const
    {
      return error;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return 1;
    }
    
    //! print some information
    void get_info () const
    {
    }

  private:
    size_type verbosity;
    const M& model;
    const S& newton;
    time_type t, dt, dtmax;
    number_type reduction;
    size_type linesearchsteps;
    Vector<number_type> u;
    Vector<number_type> unew;
    mutable bool error;
  };

  /** @brief Implementation of a general Diagonal Implicit Runge-Kutta
      method

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
      \tparam S nonlinear solver
  */
  template<class M, class S>
  class DIRK
  {
  public:

    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    /** \brief the type of a Butcher tableau */
    typedef DenseMatrix<number_type> ButcherTableau;

  private:

    //! Some default tableaus for the standard RK methods of order one
    //! to four.
    static ButcherTableau initTableau(const std::string method)
    {
      if(method.find("Implicit Euler") != std::string::npos){
        ButcherTableau butcher(2,2,0.0);
        butcher[1][1] = 1;
        butcher[0][1] = 1;
        butcher[0][0] = 1;

        return butcher;
      }
      else if(method.find("Alexander") != std::string::npos){
        ButcherTableau butcher(3,3,0.0);
        const number_type alpha = 1. - sqrt(2.)/2.;
        butcher[0][0] = alpha;
        butcher[0][1] = alpha;

        butcher[1][0] = 1.;
        butcher[1][1] = 1. - alpha;
        butcher[1][2] = alpha;

        butcher[2][1] = 1. - alpha;
        butcher[2][2] = alpha;

        return butcher;
      }
      else if(method.find("Crouzieux") != std::string::npos){
        ButcherTableau butcher(3,3,0.0);
        const number_type beta = 1./2./sqrt(3);
        butcher[0][0] = 0.5 + beta;
        butcher[0][1] = 0.5 + beta;

        butcher[1][0] = 0.5 - beta;
        butcher[1][1] = -1. / sqrt(3);
        butcher[1][2] = 0.5 + beta;

        butcher[2][1] = 0.5;
        butcher[2][2] = 0.5;

        return butcher;
      }
      else if(method.find("Midpoint Rule") != std::string::npos){
        ButcherTableau butcher(2,2,0.0);
        butcher[0][0] = 0.5;
        butcher[0][1] = 0.5;
        butcher[1][1] = 1;

        return butcher;
      }
      else if(method.find("Fractional Step Theta") != std::string::npos){
        ButcherTableau butcher(5,5,0.0);
        const number_type theta = 1 - sqrt(2.)/2.;
        const number_type alpha = 2. - sqrt(2.);
        const number_type beta = 1. - alpha;
        butcher[1][0] = theta;
        butcher[1][1] = beta * theta;
        butcher[1][2] = alpha * theta;

        butcher[2][0] = 1.-theta;
        butcher[2][1] = beta * theta;
        butcher[2][2] = alpha * (1.-theta);
        butcher[2][3] = alpha * theta;

        butcher[3][0] = 1.;
        butcher[3][1] = beta * theta;
        butcher[3][2] = alpha * (1.-theta);
        butcher[3][3] = (alpha + beta) * theta;
        butcher[3][4] = alpha * theta;

        butcher[4][1] = beta * theta;
        butcher[4][2] = alpha * (1.-theta);
        butcher[4][3] = (alpha + beta) * theta;
        butcher[4][4] = alpha * theta;

        return butcher;
      }
      else{
        HDNUM_ERROR("Order not available for Runge Kutta solver.");
      }
    }

    static int initOrder(const std::string method)
    {
      if(method.find("Implicit Euler") != std::string::npos){
        return 1;
      }
      else if(method.find("Alexander") != std::string::npos){
        return 2;
      }
      else if(method.find("Crouzieux") != std::string::npos){
        return 3;
      }
      else if(method.find("Midpoint Rule") != std::string::npos){
        return 2;
      }
      else if(method.find("Fractional Step Theta") != std::string::npos){
        return 2;
      }
      else{
        HDNUM_ERROR("Order not available for Runge Kutta solver.");
      }
    }
    

    //! class providing nonlinear problem to be solved 
    // h_n f(t_n, y_n) - y_n + y_{n-1} = 0
    class NonlinearProblem
    {
    public:
      /** \brief export size_type */
      typedef typename M::size_type size_type;
      
      /** \brief export number_type */
      typedef typename M::number_type number_type;

      //! constructor stores parameter lambda
      NonlinearProblem (const M& model_, const Vector<number_type>& yold_, 
                        typename M::time_type told_, typename M::time_type dt_,
                        const ButcherTableau & butcher_, const int rk_step_,
                        const std::vector< Vector<number_type> > & k_)
        : model(model_), yold(yold_), told(told_), 
          dt(dt_), butcher(butcher_), rk_step(rk_step_), k_old(model.size(),0)
      {
        for(int i=0; i<rk_step; ++i)
          k_old.update(butcher[rk_step][1+i] * dt, k_[i]);
      }

      //! return number of componentes for the model
      std::size_t size () const
      {
        return model.size();
      }

      //! model evaluation
      void F (const Vector<number_type>& x, Vector<number_type>& result) const
      {
        result = k_old;

        Vector<number_type> current_z(x);
        current_z.update(1.,yold);

        const number_type tnew = told + butcher[rk_step][0] * dt;

        Vector<number_type> current_k(model.size(),0.);
        model.f(tnew,current_z,current_k);
        result.update(butcher[rk_step][rk_step+1] * dt, current_k);

        result.update(-1.,x);
      }

      //! jacobian evaluation needed for implicit solvers
      void F_x (const Vector<number_type>& x, DenseMatrix<number_type>& result) const
      {
        const number_type tnew = told + butcher[rk_step][0] * dt;

        Vector<number_type> current_z(x);
        current_z.update(1.,yold);

        model.f_x(tnew,current_z,result);

        result *= dt * butcher[rk_step][rk_step+1];

        for (size_type i=0; i<model.size(); i++) result[i][i] -= number_type(1.0);
      }

      void set_told_dt (typename M::time_type told_, typename M::time_type dt_)
      {
        told = told_;
        dt = dt_;
      }

    private:
      const M& model;
      const Vector<number_type>& yold;
      typename M::time_type told;
      typename M::time_type dt;
      const ButcherTableau & butcher;
      const int rk_step;
      Vector<number_type> k_old;
    };

  public:

    //! constructor stores reference to the model and requires a
    //! butcher tableau
    DIRK (const M& model_, const S& newton_, const ButcherTableau & butcher_, const int order_)
      : verbosity(0), butcher(butcher_), model(model_), newton(newton_), 
        u(model.size()), order(order_)
    {
      model.initialize(t,u);
      dt = dtmax = 0.1; 
    }

    //! constructor stores reference to the model and sets the default
    //! butcher tableau corresponding to the given order
    DIRK (const M& model_, const S& newton_, const std::string method)
      : verbosity(0), butcher(initTableau(method)), model(model_), newton(newton_), u(model.size()),
        order(initOrder(method))
    {
      model.initialize(t,u);
      dt = dtmax = 0.1; 
    }


    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dtmax = dt_;
    }

    //! set verbosity level
    void set_verbosity (size_type verbosity_)
    {
      verbosity = verbosity_;
    }

    //! do one step
    void step ()
    {

      const size_type R = butcher.colsize()-1;

      bool reduced = false;
      error = false;
      if(verbosity>=2)
        std::cout << "DIRK: step to" << " t+dt=" << t+dt << " dt=" << dt << std::endl;

      while (1)
        {
          bool converged = true;

          // Perform R Runge-Kutta steps
          std::vector< Vector<number_type> > k;
          for(size_type i=0; i<R; ++i) {
            if (verbosity>=2) 
              std::cout << "DIRK: step nr "<< i << std::endl;

            Vector<number_type> current_z(model.size(),0.0);

            // Set starting value of k_i
            // model.f(t,u,current_k);

            // Solve nonlinear problem
            NonlinearProblem nlp(model,u,t,dt,butcher,i,k);

            newton.solve(nlp,current_z);

            converged = converged && newton.has_converged();
            if(!converged)
              break;

            current_z.update(1., u);
            const number_type t_i = t + butcher[i][0] * dt;
            Vector<number_type>current_k(model.size(),0.);
            model.f(t_i,current_z,current_k);

            k.push_back( current_k );
          }

          if (converged)
            {
              if(verbosity >= 2)
                std::cout << "DIRK finished"<< std::endl;
              
              // Update to new solution
              for(size_type i=0; i<R; ++i)
                u.update(dt*butcher[R][1+i],k[i]);

              t += dt;
              if (!reduced && dt<dtmax-1e-13)
                {
                  dt = std::min(2.0*dt,dtmax);
                  if (verbosity>0) 
                    std::cout << "DIRK: increasing time step to " << dt << std::endl;
                }
              return;
            }
          else
            {
              if (dt<1e-12)
                {
                  HDNUM_ERROR("time step too small in implicit Euler");
                  error = true;
                  break;
                }
              dt *= 0.5;
              reduced = true;
              if (verbosity>0) std::cout << "DIRK: reducing time step to " << dt << std::endl;
            }
        }
    }

    //! get current state
    bool get_error () const
    {
      return error;
    }

    //! set current state
    void set_state (time_type t_, const Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const Vector<number_type>& get_state () const
    {
      return u;
    }

    //! get current time
    time_type get_time () const
    {
      return t;
    }

    //! get dt used in last step (i.e. to compute current state)
    time_type get_dt () const
    {
      return dt;
    }

    //! return consistency order of the method
    size_type get_order () const
    {
      return order;
    }
    
    //! print some information
    void get_info () const
    {
    }

  private:
    size_type verbosity;
    const DenseMatrix<number_type> butcher; 
    const M& model;
    const S& newton;
    time_type t, dt, dtmax;
    number_type reduction;
    size_type linesearchsteps;
    Vector<number_type> u;
    int order;
    mutable bool error;
  };

  //! gnuplot output for time and state sequence
  template<class T, class N>
  inline void gnuplot (const std::string& fname, const std::vector<T> t, const std::vector<Vector<N> > u)
  {
    std::fstream f(fname.c_str(),std::ios::out);
    for (typename std::vector<T>::size_type n=0; n<t.size(); n++)
      {
        f << std::scientific << std::showpoint 
          << std::setprecision(16) << t[n];
        for (typename Vector<N>::size_type i=0; i<u[n].size(); i++)
          f << " " << std::scientific << std::showpoint 
            << std::setprecision(u[n].precision()) << u[n][i];
        f << std::endl;
      }
    f.close();
  }

  //! gnuplot output for time and state sequence
  template<class T, class N>
  inline void gnuplot (const std::string& fname, const std::vector<T> t, const std::vector<Vector<N> > u, const std::vector<T> dt)
  {
    std::fstream f(fname.c_str(),std::ios::out);
    for (typename std::vector<T>::size_type n=0; n<t.size(); n++)
      {
        f << std::scientific << std::showpoint 
          << std::setprecision(16) << t[n];
        for (typename Vector<N>::size_type i=0; i<u[n].size(); i++)
          f << " " << std::scientific << std::showpoint 
            << std::setprecision(u[n].precision()) << u[n][i];
        f << " " << std::scientific << std::showpoint 
          << std::setprecision(16) << dt[n];
        f << std::endl;
      }
    f.close();
  }

} // namespace hdnum

#endif
