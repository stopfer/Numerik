#include<math.h>

/** @brief Hodgkin Huxley model

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class HodgkinHuxley
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  //! constructor stores parameter lambda
  HodgkinHuxley () 
  : Cm(1.0),GNa(120.0), GK(36.0), Gm(0.3),
    ENa(115.0), EK(-12.0),Em(10.613)
  {}

  //! return number of componentes for the model
  std::size_t size () const
  {
    return 4;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;
    for (size_type i=0; i<3; i++) x0[i] = 0.0;
  }

  //! model evaluation
  void f (const T& t, const Vector<N>& x, Vector<N>& result) const
  {
    number_type V=x[0];
    number_type m=x[1];
    number_type h=x[2];
    number_type n=x[3];

    result[0] = (Isource(t) + GNa*m*m*m*h*(ENa-V) + GK*n*n*n*n*(EK-V) + Gm*(Em-V))/Cm;
    result[1] = alpham(V)*(1-m)-betam(V)*m;
    result[2] = alphah(V)*(1-h)-betah(V)*h;
    result[3] = alphan(V)*(1-n)-betan(V)*n;
  }

private:
  number_type Cm;
  number_type GNa, GK, Gm;
  number_type ENa, EK, Em;

  number_type alphan (number_type V) const
  {
	return (10-V)/(100.0*(exp((10-V)/10)-1));
  }

  number_type betan (number_type V) const
  {
	return 0.125*exp(-V/80);
  }

  number_type alpham (number_type V) const
  {
	return (25-V)/(10.0*(exp((25-V)/10)-1));
  }

  number_type betam (number_type V) const
  {
	return 4*exp(-V/18);
  }

  number_type alphah (number_type V) const
  {
	return 0.07*exp(-V/20);
  }

  number_type betah (number_type V) const
  {
	return 1.0/(exp((30-V)/10)+1);
  }

  number_type Isource (time_type t) const
  {
    if (t<100)
      return 10.0;
    else
      return 0.0;
  }
};
