/** @brief Explicit Euler method as an example for an ODE solver

    The ODE solver is parametrized by a model. The model also
    exports all relevant types for time and states.
    The ODE solver encapsulates the states needed for the computation.

    \tparam M the model type
*/
template<class M>
class ExplicitEuler
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

  //! constructor stores reference to the model
  ExplicitEuler (const M& model_)
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
