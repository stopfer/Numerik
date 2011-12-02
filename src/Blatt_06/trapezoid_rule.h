/*
 * trapezoid_rule.h
 *
 *  Created on: 01.12.2011
 *      Author: christopher
 */

#ifndef TRAPEZOID_RULE_H_
#define TRAPEZOID_RULE_H_

#include "hdnum.hh"
#include <gmpxx.h>
  /** @brief Implicit Trapezoid Rule

      The ODE solver is parametrized by a model. The model also
      exports all relevant types for time and states.
      The ODE solver encapsulates the states needed for the computation.

      \tparam M the model type
  */
  template<class M>
  class Trapezoid
  {
  public:
    /** \brief export size_type */
    typedef typename M::size_type size_type;

    /** \brief export time_type */
    typedef typename M::time_type time_type;

    /** \brief export number_type */
    typedef typename M::number_type number_type;

    //! constructor stores reference to the model
    Trapezoid (const M& model_)
      : model(model_), u(model.size()), f(model.size()), A(model.size(), model.size()), E(model.size(), model.size(), 0.0)
    {
      for(size_type i = 0; i < model.size(); i++) {
        E[i][i] = 1.0;
      }
    }

    //! set time step for subsequent steps
    void set_dt (time_type dt_)
    {
      dt = dt_;
    }

    //! do one step
    void step ()
    {
      mpf_set_default_prec(1024);
      // (E_n - 1/2 * h * A) * y_n = y_(n-1) + 1/2 * A * y_(n-1)
      model.f_x(t, u, A);         // get Matrix
      A *= -(1.0/2.0)*dt;         // Set Matrix for lr decomposition
      A += E;
      model.f(t, u, f);           // evaluate model
      u.update((1.0/2.0)*dt,f);   // set Vector for lr decomposition

      // lr decomposition
      hdnum::Vector<number_type> s(model.size());
      hdnum::Array<size_type> p(model.size());
      hdnum::Array<size_type> q(model.size());
      row_equilibrate(A,s);
      lr_fullpivot(A,p,q);
      apply_equilibrate(s,u);
      f = 0.0;
      permute_forward(p,u);
      solveL(A,u,u);
      solveR(A,f,u);
      permute_backward(q,f);

      u = f;            // set new state
      t += dt;          // advance time
    }

    //! set current state
    void set_state (time_type t_, const hdnum::Vector<number_type>& u_)
    {
      t = t_;
      u = u_;
    }

    //! get current state
    const hdnum::Vector<number_type>& get_state () const
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

  private:
    const M& model;
    time_type t, dt;
    hdnum::Vector<number_type> u, f;
    hdnum::DenseMatrix<number_type> A, E;
  };


#endif /* TRAPEZOID_RULE_H_ */
