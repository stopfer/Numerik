#include "nbody.hh"

/** @brief three bodies in 2d space

    \tparam T a type representing time values
    \tparam N a type representing states and f-values
*/
template<class T, class N=T>
class ThreeBody : public NBody<T,N,2>
{
public:
  /** \brief export size_type */
  typedef std::size_t size_type;

  /** \brief export time_type */
  typedef T time_type;

  /** \brief export number_type */
  typedef N number_type;

  // make a the two body model 
  ThreeBody () : NBody<T,N,2>(3)
  {
    this->m[0] = 1; 
    this->m[1] = 1; 
    this->m[2] = 1;
    this->G = 1.0;
  }

  //! set initial state including time value
  void initialize (T& t0, Vector<N>& x0) const
  {
    t0 = 0;

    number_type v(0.5);
    number_type c(1.1); // 0.99 unstable, 1.00 stable

    x0[0] = c;
    x0[1] = 0;
    x0[2] = 0;
    x0[3] = v;

    x0[4] = -0.5;
    x0[5] = 0.8660254;
    x0[6] = -0.8660254*v;
    x0[7] = -0.5*v;

    x0[8] = -0.5;
    x0[9] = -0.8660254;
    x0[10] = 0.8660254*v;
    x0[11] = -0.5*v;
    normalize(x0);
  }
};
