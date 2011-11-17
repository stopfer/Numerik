#ifndef HEUNS_METHOD_H_
#define HEUNS_METHOD_H_

/** @brief Heun's method as an example for an ODE solver

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class HeunsMethod
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

  //! constructor stores reference to the model
  HeunsMethod (const M& model_)
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

  hdnum::Vector<number_type> k1() {
    hdnum::Vector<number_type> temp(model.size());
    model.f(t,u,temp);   // evaluate model
    return temp;
  }

  hdnum::Vector<number_type> k2() {
    hdnum::Vector<number_type> temp(model.size());
    model.f(t + ((1.0/3.0) * dt),u + (k1() *= (1.0 / 3.0) * dt), temp);   // evaluate model
    return temp;
  }

  hdnum::Vector<number_type> k3() {
    hdnum::Vector<number_type> temp(model.size());
    model.f(t + ((2.0/3.0) * dt),u + (k2() *= (2.0 / 3.0) * dt), temp);   // evaluate model
    return temp;
  }

  //! do one step
  void step ()
  {
    f = (k1() += (k3() *= 3)) *= (1.0/4.0);
    u.update(dt,f);   // advance state
    t += dt;          // advance time
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
  hdnum::Vector<number_type> u;
  hdnum::Vector<number_type> f;
};


#endif /* HEUNS_METHOD_H_ */
