#include <assert.h>
#include "hdnum.hh"

using namespace hdnum;

template<class M>
class AdamsBashfort
{
public:
  /** \brief export size_type */
  typedef typename M::size_type size_type;

  /** \brief export time_type */
  typedef typename M::time_type time_type;

  /** \brief export number_type */
  typedef typename M::number_type number_type;

public:

  AdamsBashfort (const M& model_, Vector<number_type> u2_, Vector<number_type> u3_, Vector<number_type> u4_, number_type dt_)
    : model(model_), u(model.size()), u1(model.size()), u2(model.size()), u3(model.size()), u4(model.size()),
	f(model.size()), f1(model.size()), f2(model.size()), f3(model.size()), f4(model.size())
  {
    model.initialize(t,u);
    u4 = u4_;
    u3 = u3_;
    u2 = u2_;
    u1 = u;
    dt = dt_;
    tstart = t;
    t = t+3*dt;
  }

  //! do one step
  void step ()
  {
    // compute f1 to f4
    model.f(t-3*dt, u1, f1);
    f1 *= -9.0;
    model.f(t-2*dt, u2, f2);
    f2 *= 37.0;
    model.f(t-1*dt, u3, f3);
    f3 *= -59.0;
    model.f(t, u4, f4);
    f4 *= 55.0;
    // compute f
    f = f1 + f2 + f3 + f4;
    f *= dt / 24.0;

    // compute u
    u = u4 + f;

    // t and u1 to u4
    t += dt;
    u1 = u2;
    u2 = u3;
    u3 = u4;
    u4 = u;
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

private:
  const M& model;                    
  time_type t, tstart, dt;
  Vector<number_type> u;
  Vector<number_type> u1;
  Vector<number_type> u2;
  Vector<number_type> u3;
  Vector<number_type> u4;
  Vector<number_type> f;
  Vector<number_type> f1;
  Vector<number_type> f2;
  Vector<number_type> f3;
  Vector<number_type> f4;
};
